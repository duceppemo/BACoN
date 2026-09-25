"""BACoN: Bait, Assemble and Compare Nanopore reads matching a reference sequence."""

__version__ = "0.3.0"
__author__ = "Marc-Olivier Duceppe"


class BaconError(Exception):
    """A user-facing error: bad input, missing tool, or a failed external command."""
