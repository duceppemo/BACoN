"""Running external programs, with their output logged to a file."""

from __future__ import annotations

import gzip
import logging
import os
import re
import shlex
import shutil
import signal
import subprocess
import threading
from collections.abc import Sequence
from pathlib import Path

from bacon import BaconError

log = logging.getLogger(__name__)

# Conda package providing each executable, for the "missing tool" message.
PACKAGES = {
    "minimap2": "minimap2", "samtools": "samtools", "bbduk.sh": "bbmap", "filtlong": "filtlong", "flye": "flye",
    "myloasm": "myloasm", "ska": "ska2", "parsnp": "parsnp", "harvesttools": "harvesttools",
    "FastTree": "fasttree", "iqtree": "iqtree", "Bandage": "bandage",
}
# Alternative executable names, tried in order.
ALIASES = {"FastTree": ("FastTree", "fasttree"), "iqtree": ("iqtree3", "iqtree2", "iqtree")}


# Programs running now, in every thread: killed when BACoN is interrupted (see kill_running). Once stopping,
# no new program starts (a thread could otherwise start one just after the others were killed).
_RUNNING: set[subprocess.Popen] = set()
_RUNNING_LOCK = threading.Lock()
_STOPPING = threading.Event()


def kill_running() -> int:
    """Stop starting programs, and kill every program started by run() that is still running. Returns how many
    were killed. allow_programs() lifts the stop."""
    with _RUNNING_LOCK:
        _STOPPING.set()
        procs = [p for p in _RUNNING if p.poll() is None]
        _RUNNING.intersection_update(procs)  # Programs that ended (killed earlier) are forgotten
    for p in procs:
        _kill(p)
    return len(procs)


def allow_programs() -> None:
    _STOPPING.clear()


def stopping() -> bool:
    """True once kill_running() was called (until allow_programs())."""
    return _STOPPING.is_set()


def _kill(proc: subprocess.Popen) -> None:
    """Kill a program and the programs it started (each program runs in its own process group: Flye, for
    one, starts helpers that would otherwise outlive it)."""
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        proc.kill()


class ToolError(BaconError):
    """An external program failed."""


def which(name: str) -> str | None:
    for candidate in ALIASES.get(name, (name,)):
        path = shutil.which(candidate)
        if path:
            return path
    return None


def require(names: Sequence[str]) -> dict[str, str]:
    """Resolve executables, or raise one error that lists every missing one."""
    found, missing = {}, []
    for name in dict.fromkeys(names):
        path = which(name)
        if path:
            found[name] = path
        else:
            missing.append(name)
    if missing:
        packages = " ".join(sorted({PACKAGES.get(n, n) for n in missing}))
        raise BaconError(f"Required program(s) not found on PATH: {', '.join(missing)}. "
                         f"Install with: conda install -c conda-forge -c bioconda {packages}")
    return found


def _tail(path: Path, lines: int = 15) -> str:
    try:
        text = path.read_text(errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join("    " + line for line in text[-lines:])


def run(cmds: Sequence[Sequence[str]] | Sequence[str], log_file: Path, *, stdout: Path | None = None,
        cwd: Path | None = None, env: dict[str, str] | None = None, what: str = "") -> None:
    """Run a command, or a pipeline of commands, and raise ToolError if any of them fails.

    The standard error of every command, and the standard output of the last one unless it is redirected to
    `stdout` (gzipped if the name ends with .gz), are appended to `log_file`.
    """
    pipeline = [list(cmds)] if isinstance(cmds[0], str) else [list(c) for c in cmds]
    log_file.parent.mkdir(parents=True, exist_ok=True)
    full_env = {**os.environ, **env} if env else None
    command = " | ".join(shlex.join(c) for c in pipeline)
    if cwd is not None:  # The command line as it can be run again: from its working folder
        command = f"(cd {shlex.quote(str(cwd))} && {command})"
    with open(log_file, "a") as log_fh:
        log_fh.write("$ " + command + (f" > {shlex.quote(str(stdout))}" if stdout else "") + "\n")
        log_fh.flush()
        log.debug("Running: %s", command)
        procs: list[subprocess.Popen] = []
        prev = None
        try:
            for i, cmd in enumerate(pipeline):
                last = i == len(pipeline) - 1
                out = subprocess.PIPE if (not last or stdout) else log_fh
                with _RUNNING_LOCK:
                    if _STOPPING.is_set():
                        raise ToolError(f"{Path(cmd[0]).name} not started: BACoN is stopping")
                    procs.append(subprocess.Popen(cmd, stdin=prev, stdout=out, stderr=log_fh, cwd=cwd,
                                                  env=full_env, start_new_session=True))
                    _RUNNING.add(procs[-1])
                if prev is not None:
                    prev.close()  # So that the upstream process gets SIGPIPE if the downstream one dies
                prev = procs[-1].stdout
            if stdout is not None:
                gz = str(stdout).endswith(".gz")
                with (gzip.open(stdout, "wb", compresslevel=4) if gz else open(stdout, "wb")) as out_fh:
                    shutil.copyfileobj(procs[-1].stdout, out_fh, 1 << 20)
                procs[-1].stdout.close()
            codes = [p.wait() for p in procs]
        except FileNotFoundError as exc:
            for p in procs:
                _kill(p)
            raise ToolError(f"Program not found: {exc.filename}") from None
        except BaseException:
            for p in procs:
                _kill(p)
            raise
        finally:
            # A program still running stays listed: a second Ctrl-C can stop the loop above before it was killed,
            # and kill_running() must still find it.
            with _RUNNING_LOCK:
                _RUNNING.difference_update(p for p in procs if p.poll() is not None)
    failed = [(cmd, code) for cmd, code in zip(pipeline, codes) if code != 0]
    if failed:
        # In a pipeline, a program killed by SIGPIPE (-13) only reports that a later one stopped reading:
        # blame the program that failed on its own.
        cmd, code = next(((c, k) for c, k in failed if k > 0), failed[0])
        raise ToolError(f"{Path(cmd[0]).name} failed with exit code {code}{' ' + what if what else ''}; "
                        f"see {log_file}\n{_tail(log_file)}")


def version(name: str) -> str:
    """First line of a program's version output, or '' if it cannot be determined."""
    path = which(name)
    if not path:
        return ""
    args = {"FastTree": ["-expert"]}.get(name, ["--version"])
    try:
        proc = subprocess.run([path, *args], capture_output=True, text=True, timeout=60, errors="replace",
                              env={**os.environ, "QT_QPA_PLATFORM": "offscreen"})
    except (OSError, subprocess.TimeoutExpired):
        return ""
    text = proc.stdout + "\n" + proc.stderr
    if name == "FastTree":  # "Detailed usage for FastTree 2.2.0 Double precision:"
        match = re.search(r"FastTree (\d[\w.]*)", text)
        return f"FastTree {match.group(1)}" if match else ""
    for line in text.splitlines():
        line = line.strip()
        if name == "bbduk.sh" and "version" not in line.lower():
            continue
        if line and not line.startswith(("java ", "WARNING", "Warning", "/")):
            return line[:120]
    return ""
