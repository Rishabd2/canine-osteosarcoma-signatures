"""Call the R helpers in R/ with Rscript ($RSCRIPT if set, else Rscript on PATH)."""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def rscript() -> str:
    exe = os.environ.get("RSCRIPT") or shutil.which("Rscript")
    if not exe:
        raise RuntimeError("Rscript not found; set RSCRIPT or install R (see README, Reproduce)")
    return exe


def run_r(script: str, *args: str) -> str:
    cmd = [rscript(), str(ROOT / "R" / script), *map(str, args)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"{script} failed:\n{res.stdout}\n{res.stderr}")
    return res.stdout
