import gzip
from pathlib import Path

import pytest


def write_fastq(path: Path, reads: list[tuple[str, str]], gz: bool = True) -> Path:
    text = "".join(f"@{name} extra\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in reads)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    return path


@pytest.fixture
def fastq(tmp_path):
    return lambda name, reads, gz=True: write_fastq(tmp_path / name, reads, gz)


@pytest.fixture
def fasta(tmp_path):
    return lambda name, records: write_fasta(tmp_path / name, records)
