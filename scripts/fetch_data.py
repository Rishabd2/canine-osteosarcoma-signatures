"""Download every public input this replication needs and verify checksums.

Nothing fetched here is committed. Files land in ``$COS_DATA_DIR`` (default
``data/raw``). Upstream md5 sums (Zenodo, GDC) are checked when the source
publishes them; every file's sha256 is compared against ``data/checksums.sha256``
when that manifest already lists it, and written there otherwise.

Usage:
    python scripts/fetch_data.py            # fetch missing files, verify all
    python scripts/fetch_data.py --update-manifest
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import requests
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cos_sig.paths import data_dir  # noqa: E402

MANIFEST = ROOT / "data" / "checksums.sha256"

ZENODO_RECORD = "8125077"
GEO_FILES = {
    "geo/GSE238110_RawCountFile_combined.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/GSE238110/suppl/GSE238110_RawCountFile_combined.csv.gz",
    "geo/GSE238110_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/GSE238110/matrix/GSE238110_series_matrix.txt.gz",
}
_ESM = "https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-023-05208-z/MediaObjects/"
SUPP_FILES = {
    "supp/42003_2023_5208_MOESM2_ESM.pdf": _ESM + "42003_2023_5208_MOESM2_ESM.pdf",   # Supplementary Information
    "supp/42003_2023_5208_MOESM4_ESM.xlsx": _ESM + "42003_2023_5208_MOESM4_ESM.xlsx",  # Supplementary Data 1 (source data)
    "supp/42003_2023_5208_MOESM5_ESM.xlsx": _ESM + "42003_2023_5208_MOESM5_ESM.xlsx",  # Supplementary Data 2 (signatures, DEGs)
}
REF_FILES = {
    "ref/h.all.v7.5.1.symbols.gmt":
        "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/h.all.v7.5.1.symbols.gmt",
    "ref/c2.cp.reactome.v7.5.1.symbols.gmt":
        "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/c2.cp.reactome.v7.5.1.symbols.gmt",
    # NCBI CanFam4 (UU_Cfam_GSD_1.0) annotation, the reference the authors aligned to; used for gene lengths
    "ref/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gtf.gz":
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/100/685/GCF_011100685.1_UU_Cfam_GSD_1.0/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gtf.gz",
    # quanTIseq TIL10 signature matrix (Finotello et al. 2019), pinned commit
    "ref/TIL10_signature.txt":
        "https://raw.githubusercontent.com/icbi-lab/quanTIseq/{qts}/quantiseq/deconvolution/TIL10_signature.txt",
    # MCP-counter marker genes (Becht et al. 2016), pinned commit
    "ref/MCPcounter_genes.txt":
        "https://raw.githubusercontent.com/ebecht/MCPcounter/{mcp}/Signatures/genes.txt",
}
QUANTISEQ_COMMIT = "039eb16455373915ace33b3ecede09bde697d005"
MCPCOUNTER_COMMIT = "b6eac73e91c246fcff0bb1a5c68a816cd588fc48"

BIOMART = "https://www.ensembl.org/biomart/martservice"
BIOMART_QUERY = """<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">
<Dataset name="clfamiliaris_gene_ensembl" interface="default">
<Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/>
<Attribute name="hsapiens_homolog_ensembl_gene"/><Attribute name="hsapiens_homolog_associated_gene_name"/>
<Attribute name="hsapiens_homolog_orthology_type"/><Attribute name="hsapiens_homolog_orthology_confidence"/>
</Dataset></Query>"""

GDC = "https://api.gdc.cancer.gov"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def get(url: str, dest: Path, *, method: str = "GET", data=None, tries: int = 4) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(tries):
        try:
            with requests.request(method, url, data=data, stream=True, timeout=300,
                                  headers={"User-Agent": "cos-sig-replication/1.0"}) as r:
                r.raise_for_status()
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
            tmp.replace(dest)
            return
        except requests.RequestException as exc:  # pragma: no cover - network
            if attempt == tries - 1:
                raise
            print(f"  retry {attempt + 1} for {url}: {exc}")
            time.sleep(5 * (attempt + 1))


def fetch_zenodo(base: Path) -> dict[str, str]:
    rec = requests.get(f"https://zenodo.org/api/records/{ZENODO_RECORD}", timeout=60).json()
    md5s = {}
    for f in rec["files"]:
        rel = f"zenodo/{f['key']}"
        dest = base / rel
        expect = f["checksum"].split(":", 1)[1]
        if not dest.exists() or md5(dest) != expect:
            print(f"fetch {rel}")
            get(f["links"]["self"], dest)
        if md5(dest) != expect:
            raise SystemExit(f"md5 mismatch for {rel}")
        md5s[rel] = expect
    return md5s


def fetch_gdc(base: Path) -> None:
    filt = {"op": "and", "content": [
        {"op": "=", "content": {"field": "cases.project.project_id", "value": "TARGET-OS"}},
        {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
        {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
        {"op": "=", "content": {"field": "access", "value": "open"}}]}
    fields = "file_id,file_name,md5sum,file_size,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type"
    r = requests.post(f"{GDC}/files", json={"filters": filt, "size": 2000, "fields": fields}, timeout=120).json()
    hits = r["data"]["hits"]
    rows = []
    for h in hits:
        case = h["cases"][0]
        rows.append({"file_id": h["file_id"], "file_name": h["file_name"], "md5sum": h["md5sum"],
                     "case": case["submitter_id"], "sample": case["samples"][0]["submitter_id"],
                     "sample_type": case["samples"][0]["sample_type"]})
    (base / "gdc").mkdir(parents=True, exist_ok=True)
    (base / "gdc" / "file_manifest.json").write_text(json.dumps(sorted(rows, key=lambda x: x["file_id"]), indent=1))
    print(f"GDC: {len(rows)} TARGET-OS STAR-Counts files")
    for row in rows:
        dest = base / "gdc" / "star_counts" / f"{row['file_id']}.tsv"
        if dest.exists() and md5(dest) == row["md5sum"]:
            continue
        get(f"{GDC}/data/{row['file_id']}", dest)
        if md5(dest) != row["md5sum"]:
            raise SystemExit(f"md5 mismatch for GDC file {row['file_id']}")
    # Clinical: GDC case records for the same project (demographics, diagnoses, follow-ups)
    cfields = ",".join([
        "submitter_id", "demographic.vital_status", "demographic.days_to_death",
        "diagnoses.days_to_last_follow_up", "diagnoses.age_at_diagnosis",
        "diagnoses.tumor_stage", "diagnoses.ajcc_clinical_m", "diagnoses.metastasis_at_diagnosis",
        "diagnoses.progression_or_recurrence", "diagnoses.days_to_recurrence",
        "follow_ups.days_to_follow_up", "follow_ups.progression_or_recurrence",
    ])
    cf = {"op": "=", "content": {"field": "project.project_id", "value": "TARGET-OS"}}
    dest = base / "gdc" / "clinical_cases.json"
    if not dest.exists():
        cr = requests.post(f"{GDC}/cases", json={"filters": cf, "size": 5000, "fields": cfields}, timeout=120).json()
        dest.write_text(json.dumps(cr["data"]["hits"], indent=1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--update-manifest", action="store_true",
                    help="rewrite data/checksums.sha256 from the files on disk")
    ap.add_argument("--skip-gdc", action="store_true")
    args = ap.parse_args()
    base = data_dir()
    base.mkdir(parents=True, exist_ok=True)
    print(f"data dir: {base}")

    fetch_zenodo(base)
    subs = {"qts": QUANTISEQ_COMMIT, "mcp": MCPCOUNTER_COMMIT}
    for rel, url in {**GEO_FILES, **SUPP_FILES, **REF_FILES}.items():
        dest = base / rel
        if not dest.exists():
            print(f"fetch {rel}")
            get(url.format(**subs), dest)
    dest = base / "ref" / "ensembl_dog_human_orthologs.tsv"
    if not dest.exists():
        print("fetch Ensembl BioMart dog-human orthologs")
        get(BIOMART + "?query=" + quote(BIOMART_QUERY), dest)
        info = requests.get("https://rest.ensembl.org/info/data?content-type=application/json", timeout=60).json()
        (base / "ref" / "ensembl_release.json").write_text(json.dumps(info))
    if not args.skip_gdc:
        fetch_gdc(base)

    # sha256 manifest over everything except the per-file GDC downloads (those are md5-checked above)
    tracked = sorted(p for p in base.rglob("*") if p.is_file() and "star_counts" not in p.parts
                     and not p.name.startswith("._") and not p.name.endswith(".part"))
    current = {str(p.relative_to(base)): sha256(p) for p in tracked}
    known = {}
    if MANIFEST.exists():
        for line in MANIFEST.read_text().splitlines():
            if line.strip():
                digest, rel = line.split(None, 1)
                known[rel.strip()] = digest
    bad = [rel for rel, d in current.items() if rel in known and known[rel] != d]
    if bad and not args.update_manifest:
        print("sha256 differs from data/checksums.sha256 for:")
        for rel in bad:
            print("  ", rel)
        print("Upstream content changed. Rerun with --update-manifest if that is expected.")
        raise SystemExit(1)
    if args.update_manifest or not MANIFEST.exists() or set(current) - set(known):
        merged = {**known, **current} if not args.update_manifest else current
        MANIFEST.write_text("".join(f"{d}  {rel}\n" for rel, d in sorted(merged.items())))
        print(f"wrote {MANIFEST.relative_to(ROOT)} ({len(merged)} files)")
    print(f"verified {len(current)} files")


if __name__ == "__main__":
    main()
