"""Published signatures from Mannheimer et al. 2023 and helpers to freeze derived ones."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

# Supplementary Table 3 (Supplementary Information PDF, page 3), transcribed verbatim.
# The main text gives sizes ISA1 30, ISA2 19, ISA3 39, ISA4 15; the table lists 29, 19, 15, 39,
# i.e. the ISA3/ISA4 sizes in the text are swapped relative to the table and to Fig. 7c,d.
PUBLISHED_ISA = {
    "ISA-1": "A4GALT ABCA8 ABI3BP ACKR2 APOA1 ARNTL2 CFD CFH CHRDL1 CXCL14 GPRC5A LSP1 MEDAG MEOX1 MFAP4 "
             "MFAP5 MYH11 NRCAM PI15 PI16 PRG4 RARRES2 SELP SEMA3C SEMA3G SLC14A1 SLC46A2 THBS4 TMEM176A".split(),
    "ISA-2": "ACTA1 CKM DES KLHL41 MYBPC1 MYH1 MYH2 MYL1 MYL2 NEB NRAP PDK4 PDLIM3 PYGM RYR1 TMOD4 TNNT1 "
             "TNNT3 TXLNB".split(),
    "ISA-3": "CD2 CD3E CD6 CLEC2D CXCR6 GBP1 GBP5 GIMAP7 GIMAP8 GZMB GZMK IL2RB IL2RG SECTM1 TBC1D10C".split(),
    "ISA-4": "APCDD1L AQP5 BARX1 BPIFB4 CACNG7 COL11A2 COL2A1 DIPK1C DLX1 DPYSL5 DUOX1 FADS6 FAT3 FGFR3 "
             "GRIA2 GRM4 HORMAD1 HOXB13 HOXC13 IGF2BP1 IHH IRX1 LHX2 LINGO1 NDNF NUP62CL PLPPR2 PPP1R1B "
             "RAB27B SCML2 SLC13A5 SLC36A2 SLC8A3 SLC9A2 SMPD3 SYT13 TBC1D16 UNC5C ZIC5".split(),
}


def published_gs(xlsx: Path) -> dict[str, list[str]]:
    """GS-1 and GS-2 from Supplementary Data 2 (sheet 'Gene Signature').

    The sheet has GS-1 (25), GS-2 (156) and 'DOG Only' (DLA-DQA1, DLA-DQB1). The
    paper says 27 GS-1 genes were present in DOG2 and that two canine-specific genes
    were absent from TARGET; Fig. 2a lists DLA-DQA1/DLA-DQB1 among the GS-1 rows.
    So the canine GS-1 is the 25 plus the two DLA genes, and the human GS-1 is the 25.
    """
    df = pd.read_excel(xlsx, sheet_name="Gene Signature")
    gs1 = df["GS-1"].dropna().astype(str).tolist()
    dog_only = df["DOG Only"].dropna().astype(str).tolist()
    gs2 = df["GS-2"].dropna().astype(str).tolist()
    return {"GS-1": gs1 + dog_only, "GS-1-human": gs1, "GS-2": gs2}


def freeze(signatures: dict[str, list[str]], path: Path, provenance: dict) -> str:
    """Write signatures to JSON and return the sha256 of the canonical gene content."""
    canon = json.dumps({k: sorted(v) for k, v in sorted(signatures.items())}, sort_keys=True)
    digest = hashlib.sha256(canon.encode()).hexdigest()
    path.write_text(json.dumps({"signatures": {k: sorted(v) for k, v in signatures.items()},
                                "sha256": digest, "provenance": provenance}, indent=1))
    return digest


def load_frozen(path: Path) -> dict[str, list[str]]:
    """Load frozen signatures and verify the stored hash (raises on any edit)."""
    blob = json.loads(Path(path).read_text())
    canon = json.dumps({k: sorted(v) for k, v in sorted(blob["signatures"].items())}, sort_keys=True)
    if hashlib.sha256(canon.encode()).hexdigest() != blob["sha256"]:
        raise RuntimeError(f"frozen signatures in {path} do not match their recorded sha256")
    return blob["signatures"]
