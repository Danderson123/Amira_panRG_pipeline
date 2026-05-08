#!/usr/bin/env python3
"""
Run Jellyfish to count canonical 31-mers for each FASTA in an assembly directory.

Writes a CSV with:
  column 1: fasta filename
  column 2: number of DISTINCT 31-mers with count > 1

Requirements:
- jellyfish available on PATH
- FASTA files (optionally gzipped) in ASSEMBLY_DIR

Notes:
- Uses canonical k-mers via jellyfish `-C`.
- Counts DISTINCT k-mers with count > 1 (i.e., number of lines in `jellyfish dump -c`
  whose count > 1), not the total occurrences.
"""

from __future__ import annotations

import csv
import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, List
from tqdm import tqdm

# =========================
# HARD-CODED CONFIG
# =========================
ASSEMBLY_DIR = "/hps/nobackup/iqbal/dander/Amira_panRG_pipeline_test/E_coli_kmer_vs_genemer"
OUTPUT_CSV = "jellyfish_31mer_count_gt1.csv"

K = 31
THREADS = 16

# Jellyfish hash size (-s). If too small, jellyfish can fail.
# You can set this to something larger (e.g., "500M", "2G") depending on genome size/complexity.
HASH_SIZE = "200M"

# FASTA filename patterns to include
FASTA_GLOBS = ["*.fa", "*.fna", "*.fasta", "*.faa", "*.fa.gz", "*.fna.gz", "*.fasta.gz", "*.faa.gz"]

# If True, uses a streaming named pipe for gzipped FASTAs to avoid writing full decompressed files.
STREAM_GZ_VIA_FIFO = True
# =========================


def run(cmd: List[str], *, check: bool = True, capture: bool = False, text: bool = True):
    """Small subprocess helper."""
    if capture:
        return subprocess.run(cmd, check=check, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=text)
    return subprocess.run(cmd, check=check)


def is_gz(path: Path) -> bool:
    return path.suffix.lower() == ".gz"


def iter_fastas(dir_path: Path) -> List[Path]:
    files: List[Path] = []
    for pat in FASTA_GLOBS:
        files.extend(dir_path.glob(pat))
    # de-dup and sort
    return sorted(set(p for p in files if p.is_file()))


def count_distinct_kmers_gt1(jf_path: Path) -> int:
    """
    Parse `jellyfish dump -c` output and count how many distinct kmers have count > 1.
    Output format is typically: "<kmer> <count>" per line.
    """
    # Use streaming dump to avoid large intermediate files.
    proc = subprocess.Popen(
        ["jellyfish", "dump", "-c", str(jf_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
        universal_newlines=True,
    )

    n = 0
    assert proc.stdout is not None
    for line in proc.stdout:
        line = line.strip()
        if not line:
            continue
        # split on whitespace
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            c = int(parts[-1])
        except ValueError:
            continue
        if c > 1:
            n += 1

    _, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"`jellyfish dump` failed:\n{err}")
    return n


def jellyfish_count_for_fasta(fasta_path: Path, workdir: Path) -> int:
    """
    Run jellyfish count (canonical 31-mers) for one FASTA and return distinct 31-mers with count>1.
    """
    jf_out = workdir / (fasta_path.name.replace(os.sep, "_") + f".k{K}.jf")

    if is_gz(fasta_path) and STREAM_GZ_VIA_FIFO:
        # Create a FIFO and stream decompression into it
        fifo = workdir / (fasta_path.stem + ".fifo")  # .stem removes .gz
        try:
            os.mkfifo(fifo)

            # Start decompressor writing to fifo
            # (python gzip -> fifo) keeps it portable without shell pipelines
            import threading

            def writer():
                with gzip.open(fasta_path, "rb") as fin, open(fifo, "wb", buffering=0) as fout:
                    shutil.copyfileobj(fin, fout)

            t = threading.Thread(target=writer, daemon=True)
            t.start()

            # jellyfish reads from fifo
            run([
                "jellyfish", "count",
                "-m", str(K),
                "-s", str(HASH_SIZE),
                "-t", str(THREADS),
                "-C",
                "-o", str(jf_out),
                str(fifo),
            ])

        finally:
            # Best-effort cleanup
            try:
                if fifo.exists():
                    fifo.unlink()
            except Exception:
                pass
    else:
        # Non-gz OR we choose not to stream: jellyfish reads the file directly
        run([
            "jellyfish", "count",
            "-m", str(K),
            "-s", str(HASH_SIZE),
            "-t", str(THREADS),
            "-C",
            "-o", str(jf_out),
            str(fasta_path),
        ])

    return count_distinct_kmers_gt1(jf_out)


def main() -> None:
    # Basic sanity check
    try:
        run(["jellyfish", "--version"], capture=True)
    except Exception as e:
        raise SystemExit(f"ERROR: jellyfish not found or not runnable on PATH: {e}")

    assembly_dir = Path(ASSEMBLY_DIR)
    if not assembly_dir.is_dir():
        raise SystemExit(f"ERROR: ASSEMBLY_DIR is not a directory: {assembly_dir}")

    fastas = iter_fastas(assembly_dir)
    if not fastas:
        raise SystemExit(f"ERROR: no FASTA files found in {assembly_dir} with patterns: {FASTA_GLOBS}")

    rows = []
    with tempfile.TemporaryDirectory(prefix="jellyfish_work_") as td:
        workdir = Path(td)
        for fasta in tqdm(fastas):
            try:
                n_gt1 = jellyfish_count_for_fasta(fasta, workdir)
            except subprocess.CalledProcessError as e:
                raise SystemExit(f"ERROR: jellyfish failed on {fasta}:\n{e}")
            rows.append((fasta.name, n_gt1))

    # Write CSV
    out_path = Path(OUTPUT_CSV)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["filename", "distinct_31mers_count_gt1"])
        w.writerows(rows)

    # Optional: print a small summary
    print(f"Wrote: {out_path.resolve()}")
    print(f"Processed {len(rows)} FASTA file(s).")


if __name__ == "__main__":
    main()
