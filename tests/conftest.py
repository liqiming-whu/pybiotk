import os
import pytest
from collections import namedtuple

from pathlib import Path
import subprocess
import sys

data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

ROOT = Path(__file__).parent.parent
SRC = ROOT / "src"
BX_DIR = SRC / "pybiotk" / "bx"


def _extension_exists(pattern):
    return any(BX_DIR.glob(pattern))


def _build_needed():
    patterns = ["bitset*.so", "cluster*.so", "intersection*.so"]

    if not all(_extension_exists(pattern) for pattern in patterns):
        return True

    extension_files = [
        path
        for pattern in patterns
        for path in BX_DIR.glob(pattern)
    ]

    source_files = [
        *ROOT.glob("src/bx/*.pyx"),
        *ROOT.glob("src/bx/*.c"),
        ROOT / "setup.py",
        ROOT / "pyproject.toml",
    ]

    newest_source = max(path.stat().st_mtime for path in source_files)
    oldest_extension = min(path.stat().st_mtime for path in extension_files)

    return newest_source > oldest_extension


def pytest_sessionstart(session):
    if not _build_needed():
        return

    BX_DIR.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [sys.executable, "setup.py", "build_ext", "--inplace"],
        cwd=ROOT,
        check=True,
    )
