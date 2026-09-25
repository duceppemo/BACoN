"""Backward-compatible entry point: `python bacon.py ...` is the same as `bacon ...`."""

import sys

from bacon.cli import main

if __name__ == "__main__":
    sys.exit(main())
