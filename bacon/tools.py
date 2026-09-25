"""Running external programs, with their output logged to a file."""

from __future__ import annotations

import gzip
import logging
import os
import re
import shlex
import shutil
import subprocess
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
    with open(log_file, "a") as log_fh:
        log_fh.write("$ " + " | ".join(shlex.join(c) for c in pipeline)
                     + (f" > {stdout}" if stdout else "") + "\n")
        log_fh.flush()
        log.debug("Running: %s", " | ".join(shlex.join(c) for c in pipeline))
        procs: list[subprocess.Popen] = []
        prev = None
        try:
            for i, cmd in enumerate(pipeline):
                last = i == len(pipeline) - 1
                out = subprocess.PIPE if (not last or stdout) else log_fh
                procs.append(subprocess.Popen(cmd, stdin=prev, stdout=out, stderr=log_fh, cwd=cwd, env=full_env))
                if prev is not None:
                    prev.close()  # So that the upstream process gets SIGPIPE if the downstream one dies
                prev = procs[-1].stdout
            if stdout is not None:
                opener = gzip.open if str(stdout).endswith(".gz") else open
                with opener(stdout, "wb") as out_fh:
                    shutil.copyfileobj(procs[-1].stdout, out_fh, 1 << 20)
                procs[-1].stdout.close()
            codes = [p.wait() for p in procs]
        except FileNotFoundError as exc:
            for p in procs:
                p.kill()
            raise ToolError(f"Program not found: {exc.filename}") from None
        except BaseException:
            for p in procs:
                p.kill()
            raise
    for cmd, code in zip(pipeline, codes):
        if code != 0:
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
