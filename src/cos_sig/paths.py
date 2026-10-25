"""Filesystem locations. Large inputs live outside git, under COS_DATA_DIR."""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def data_dir() -> Path:
    """Raw-input directory: $COS_DATA_DIR if set, else <repo>/data/raw."""
    env = os.environ.get("COS_DATA_DIR")
    return Path(env).expanduser().resolve() if env else ROOT / "data" / "raw"


def work_dir() -> Path:
    """Intermediate (derived, uncommitted) matrices: $COS_WORK_DIR or <repo>/data/work."""
    env = os.environ.get("COS_WORK_DIR")
    p = Path(env).expanduser().resolve() if env else ROOT / "data" / "work"
    p.mkdir(parents=True, exist_ok=True)
    return p


RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
