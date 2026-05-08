#!/usr/bin/env python3
"""
Aggregate CDS locus tags (with orientation) per contig for a directory of GFF/GFF3 files.

Hardcoded configuration is at the top of the file.

Output structure:
{
  "file1.gff": {
      "contigA": ["+LOCUS_0001", "-LOCUS_0002", ...],   # ordered by genomic position
      "contigB": [...],
  },
  "file2.gff3": { ... }
}

Notes:
- Only CDS features are collected.
- Orientation is taken from the GFF strand column: '+' or '-'. If missing/unknown, '?' is used.
- Locus tag is taken from attributes in this priority order:
    locus_tag, locus, ID, Name, Parent
  (If none found, a fallback "CDS_<n>" is used per file.)
- Ordering is by (start, end, locus_tag) within each contig.
"""

from __future__ import annotations

import glob
import json
import os
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
from tqdm import tqdm

# =========================
# HARD CODED CONFIGURATION
# =========================
INPUT_DIR = "/path/to/gff_dir"   # directory containing GFF/GFF3 files
GLOB_PATTERN = "*.gff*"          # glob pattern relative to INPUT_DIR
WRITE_JSON = True                # write output JSON file
OUTPUT_JSON = "aggregated_gffs.json"
PRINT_TO_SCREEN = True           # print dictionary to stdout


# =========================
# PARSING HELPERS
# =========================
def parse_gff_attributes(attr_field: str) -> Dict[str, str]:
    """
    Parse GFF3 attributes field into a dict.
    Supports typical 'key=value;key2=value2' format.
    Also tolerates older 'key value' / 'key:value' styles in a best-effort way.
    """
    attrs: Dict[str, str] = {}
    s = (attr_field or "").strip()
    if not s or s == ".":
        return attrs

    parts = [p for p in s.split(";") if p]
    for p in parts:
        p = p.strip()
        if not p:
            continue

        if "=" in p:
            k, v = p.split("=", 1)
        elif ":" in p:
            k, v = p.split(":", 1)
        elif " " in p:
            k, v = p.split(" ", 1)
        else:
            k, v = p, ""

        k = k.strip()
        v = v.strip()

        # if multiple values, take first
        if "," in v:
            v = v.split(",", 1)[0].strip()

        if k:
            attrs[k] = v

    return attrs


def pick_locus_tag(attrs: Dict[str, str]) -> Optional[str]:
    """
    Choose a locus tag-like identifier from attributes.
    Priority:
      locus_tag, locus, ID, Name, Parent
    """
    for key in ("locus_tag", "locus", "ID", "Name", "Parent"):
        val = attrs.get(key)
        if val:
            return val
    return None


# =========================
# CORE LOGIC
# =========================
def aggregate_gff_file(path: str) -> Dict[str, List[str]]:
    """
    For a single GFF file, return:
      { contig: [<strand><locus_tag>, ...] } ordered by genomic coordinates.
    """
    records: List[Tuple[str, int, int, str, str]] = []
    fallback_counter = 0

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) < 9:
                continue

            seqid, source, ftype, start_s, end_s, score, strand, phase, attrs_s = cols[:9]

            if ftype != "CDS":
                continue

            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue

            strand = strand.strip()
            if strand not in ("+", "-"):
                strand = "?"

            attrs = parse_gff_attributes(attrs_s)
            locus = pick_locus_tag(attrs)
            if not locus:
                fallback_counter += 1
                locus = f"CDS_{fallback_counter}"

            records.append((seqid, start, end, strand, locus))

    by_contig: Dict[str, List[Tuple[int, int, str, str]]] = defaultdict(list)
    for contig, start, end, strand, locus in records:
        by_contig[contig].append((start, end, strand, locus))

    result: Dict[str, List[str]] = {}
    for contig, feats in by_contig.items():
        feats.sort(key=lambda x: (x[0], x[1], x[3]))
        result[contig] = [f"{strand}{locus}" for (start, end, strand, locus) in feats]

    return result


def aggregate_gff_dir(input_dir: str, pattern: str) -> Dict[str, Dict[str, List[str]]]:
    """
    Glob input_dir/pattern and aggregate all matching GFF/GFF3 files into:
      { filename: { contig: [<strand><locus_tag>, ...], ... }, ... }
    """
    glob_path = "/hps/nobackup/iqbal/dander/Amira_panRG_pipeline_test/Escherichia_coli_panRG_kmer_vs_genemer/bakta_annotated_assemblies/*/*.gff3"
    paths = sorted(glob.glob(glob_path))

    aggregated: Dict[str, Dict[str, List[str]]] = {}
    for p in tqdm(paths):
        if not os.path.isfile(p):
            continue
        key = os.path.basename(p)
        aggregated[key] = aggregate_gff_file(p)

    return aggregated


# =========================
# MAIN EXECUTION
# =========================
def main() -> None:
    aggregated = aggregate_gff_dir(INPUT_DIR, GLOB_PATTERN)

    if PRINT_TO_SCREEN:
        print(json.dumps(aggregated, indent=2, sort_keys=True))

    if WRITE_JSON:
        with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
            json.dump(aggregated, f, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
