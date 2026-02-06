#!/usr/bin/env python3
"""
extract_v4.py
-------------
Extract V4 region from paired-end MiSeq data by primer-based trimming.

R1 logic:
    [... anything ...] [515F] [V4 ─── keep this ───>]
    → search for 515F anywhere in R1, remove everything up to and including it.

R2 logic:
    [806R] [V4 ─── keep this ───]
    → trim 806R from 5' end if found
    → if 806R not found, assume primers were already removed and use R2 as-is
    → truncate to 150bp

Default primers:
    515F          GTGYCAGCMGCCGCGGTAA      (19 bp)
    806R          GACTACNVGGGTWTCTAAT      (19 bp)

Usage:
    python extract_v4.py SAMPLE_NAME out_dir
    python extract_v4.py SAMPLE_NAME out_dir --mismatches 3

    Input  : ./SAMPLE_NAME_R1.fastq.gz  ./SAMPLE_NAME_R2.fastq.gz
    Output : out_dir/SAMPLE_NAME_trimmed_R1.fastq.gz  ...trimmed_R2.fastq.gz
"""

import gzip
import argparse
import os
import sys
import time
import logging
from typing import Optional

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─── Defaults ────────────────────────────────────────────────────────────────
V4_FW      = "GTGYCAGCMGCCGCGGTAA"       # 515F (19 bp)
V4_RV      = "GACTACNVGGGTWTCTAAT"       # 806R (19 bp)
R2_TRUNC_LEN = 150                        # Truncate R2 to this length

# ─── IUPAC ───────────────────────────────────────────────────────────────────
IUPAC = {
    "A":{"A"}, "C":{"C"}, "G":{"G"}, "T":{"T"},
    "R":{"A","G"}, "Y":{"C","T"}, "S":{"G","C"}, "W":{"A","T"},
    "K":{"G","T"}, "M":{"A","C"},
    "B":{"C","G","T"}, "D":{"A","G","T"}, "H":{"A","C","T"}, "V":{"A","C","G"},
    "N":{"A","C","G","T"},
}


# ─── Primer search ───────────────────────────────────────────────────────────
def find_primer(seq: str, primer: str, max_mm: int, search_range: Optional[int] = None) -> int:
    """
    Slide primer over seq[0 : search_range] and return the START index of the
    best match if it has <= max_mm mismatches.  Returns -1 if not found.

    search_range=None  → scan the whole read.
    """
    p_len = len(primer)
    limit = (len(seq) - p_len + 1) if search_range is None else min(search_range, len(seq) - p_len + 1)

    best_pos = -1
    best_mm  = max_mm + 1          # lower is better

    for start in range(limit):
        mm = 0
        for i in range(p_len):
            if seq[start + i].upper() not in IUPAC.get(primer[i].upper(), set()):
                mm += 1
                if mm >= best_mm:   # already worse than current best → skip
                    break
        if mm < best_mm:
            best_mm  = mm
            best_pos = start
            if mm == 0:             # perfect match → no need to keep scanning
                break

    return best_pos if best_mm <= max_mm else -1


# ─── FASTQ I/O ───────────────────────────────────────────────────────────────
def open_fq(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def iter_fastq(fh):
    while True:
        h = fh.readline().strip()
        if not h:
            break
        s = fh.readline().strip()
        _  = fh.readline()
        q = fh.readline().strip()
        yield h, s, q


# ─── Main ────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Extract V4 region from paired-end FASTQ by primer search."
    )
    parser.add_argument("sample",          help="Sample name (input: sample_R1.fastq.gz / sample_R2.fastq.gz)")
    parser.add_argument("outdir",          help="Output directory")
    parser.add_argument("--fw",            default=V4_FW,
                        help=f"V4 forward primer 515F (default: {V4_FW})")
    parser.add_argument("--rv",            default=V4_RV,
                        help=f"V4 reverse primer 806R (default: {V4_RV})")
    parser.add_argument("--mismatches",    type=int, default=2,
                        help="Max mismatches allowed per primer (default: 2)")
    parser.add_argument("--r2-rv-window",  type=int, default=30,
                        help="Search window for 806R on R2 5' end (default: 30 bp)")
    parser.add_argument("--r2-trunc-len",  type=int, default=R2_TRUNC_LEN,
                        help=f"Truncate R2 to this length after primer removal (default: {R2_TRUNC_LEN} bp)")

    args = parser.parse_args()

    # ── derive paths from sample name ──
    r1_path = f"{args.sample}_R1.fastq.gz"
    r2_path = f"{args.sample}_R2.fastq.gz"

    for p in (r1_path, r2_path):
        if not os.path.isfile(p):
            log.error(f"File not found: {p}")
            sys.exit(1)
    os.makedirs(args.outdir, exist_ok=True)

    fw      = args.fw
    rv      = args.rv
    r2_trunc = args.r2_trunc_len

    log.info(f"Sample            : {args.sample}")
    log.info(f"V4 FW primer      : {fw}  ({len(fw)} bp)")
    log.info(f"V4 RV primer      : {rv}  ({len(rv)} bp)")
    log.info(f"Max mismatches    : {args.mismatches}")
    log.info(f"R2 806R window    : {args.r2_rv_window} bp")
    log.info(f"R2 truncate len   : {r2_trunc} bp")
    log.info(f"Input             : {r1_path}  /  {r2_path}")

    out_r1 = os.path.join(args.outdir, f"{args.sample}_trimmed_R1.fastq.gz")
    out_r2 = os.path.join(args.outdir, f"{args.sample}_trimmed_R2.fastq.gz")

    # ── counters ──
    total = kept = 0
    r1_no_fw = r2_no_rv = r1_empty = r2_empty = 0
    t0 = time.time()

    with open_fq(r1_path) as fh1, open_fq(r2_path) as fh2, \
         gzip.open(out_r1, "wt") as oh1, gzip.open(out_r2, "wt") as oh2:

        for (h1, s1, q1), (h2, s2, q2) in zip(iter_fastq(fh1), iter_fastq(fh2)):
            total += 1

            # ── R1: find 515F anywhere, keep everything after it ──
            pos_fw = find_primer(s1, fw, args.mismatches)
            if pos_fw < 0:
                r1_no_fw += 1
                continue
            trim_start_r1 = pos_fw + len(fw)
            s1_out = s1[trim_start_r1:]
            q1_out = q1[trim_start_r1:]
            if not s1_out:
                r1_empty += 1
                continue

            # ── R2: find 806R near 5' end, trim it, then truncate to r2_trunc bp ──
            # If primer not found, assume primers were already removed and use R2 as-is
            pos_rv = find_primer(s2, rv, args.mismatches, search_range=args.r2_rv_window)
            if pos_rv < 0:
                # Primer not found - assume already trimmed, use R2 as-is
                r2_no_rv += 1
                s2_trimmed = s2
                q2_trimmed = q2
            else:
                trim_start_r2 = pos_rv + len(rv)
                s2_trimmed = s2[trim_start_r2:]
                q2_trimmed = q2[trim_start_r2:]

            # ── R2: truncate to specified length ──
            s2_out = s2_trimmed[:r2_trunc]
            q2_out = q2_trimmed[:r2_trunc]

            if not s2_out:
                r2_empty += 1
                continue

            # ── write ──
            oh1.write(f"{h1}\n{s1_out}\n+\n{q1_out}\n")
            oh2.write(f"{h2}\n{s2_out}\n+\n{q2_out}\n")
            kept += 1

            if total % 100_000 == 0:
                log.info(f"  ... {total:,} pairs, {kept:,} kept")

    elapsed = time.time() - t0

    log.info("─" * 58)
    log.info(f"  Total pairs processed       : {total:>10,}")
    log.info(f"  R1 — 515F not found         : {r1_no_fw:>10,}  ({r1_no_fw/total*100:5.2f}%)")
    log.info(f"  R1 — empty after trim       : {r1_empty:>10,}  ({r1_empty/total*100:5.2f}%)")
    log.info(f"  R2 — 806R not found (5')    : {r2_no_rv:>10,}  ({r2_no_rv/total*100:5.2f}%)  (used as-is)")
    log.info(f"  R2 — empty after trim       : {r2_empty:>10,}  ({r2_empty/total*100:5.2f}%)")
    log.info(f"  R2 — truncated to           : {r2_trunc:>10,} bp")
    log.info(f"  Pairs kept                  : {kept:>10,}  ({kept/total*100:5.2f}%)")
    log.info(f"  Output                      : {out_r1}")
    log.info(f"                                {out_r2}")
    log.info(f"  Elapsed                     : {elapsed:>9.1f} s")
    log.info("─" * 58)


if __name__ == "__main__":
    main()
