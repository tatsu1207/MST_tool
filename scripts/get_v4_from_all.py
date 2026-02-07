#!/usr/bin/env python3
"""
get_v4_from_all.py
------------------
Extract V4 region from paired-end MiSeq data by primer-based trimming.
Supports multiple amplicon types: v34, v4, v45.

Amplicon configurations:
    v34 (V3-V4): R1 has 515F anywhere (after V3), R2 has 806R at 5' end
    v4  (V4):    R1 has 515F at 5' end, R2 has 806R at 5' end
    v45 (V4-V5): R1 has 515F at 5' end, R2 has 926R at 5' end

R1 logic:
    v34: [... V3 ...] [515F] [V4 ─── keep this ───>]
         → search for 515F anywhere in R1, remove everything up to and including it.
    v4:  [515F] [V4 ─── keep this ───>]
         → search for 515F at 5' end only (first 30bp), remove it.
    v45: [515F] [V4 ─── keep ───] [806R] [V5 ─── discard ───>]
         → remove 515F, then find 806R and truncate before it (removes V5).

R2 logic:
    v34/v4: [806R] [V4 ─── keep this ───] → truncate to 150bp
    v45:    [926R] [V5 ─ discard ─] [806R_rc] [V4 ─── keep ───>]
            → remove 926R, find 806R_revcomp, keep from there (skips V5).

Primers:
    515F          GTGYCAGCMGCCGCGGTAA      (19 bp) - V4 forward
    806R          GACTACNVGGGTWTCTAAT      (19 bp) - V4 reverse
    926R          CCGYCAATTYMTTTRAGTTT     (20 bp) - V5 reverse

Usage:
    python get_v4_from_all.py SAMPLE_NAME out_dir
    python get_v4_from_all.py SAMPLE_NAME out_dir --amplicon v4
    python get_v4_from_all.py SAMPLE_NAME out_dir --amplicon v45 --primers trimmed
    python get_v4_from_all.py SAMPLE_NAME out_dir --primers auto --mismatches 3

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

# ─── Primer sequences ─────────────────────────────────────────────────────────
PRIMER_515F = "GTGYCAGCMGCCGCGGTAA"       # 515F (19 bp) - V4 forward
PRIMER_806R = "GACTACNVGGGTWTCTAAT"       # 806R (19 bp) - V4 reverse
PRIMER_806R_RC = "ATTAGAWACCCBNGTAGTC"    # 806R reverse complement (for R2 in V4-V5)
PRIMER_926R = "CCGYCAATTYMTTTRAGTTT"      # 926R (20 bp) - V5 reverse

# ─── Amplicon configurations ──────────────────────────────────────────────────
AMPLICON_CONFIG = {
    "v34": {
        "description": "V3-V4 amplicons (515F anywhere in R1, 806R at R2 5' end)",
        "fw_primer": PRIMER_515F,
        "rv_primer": PRIMER_806R,
        "r1_search_range": None,   # Search anywhere in R1
        "r2_trunc_len": 150,
        "extract_v4_only": False,  # Already V4 only
    },
    "v4": {
        "description": "V4 amplicons (515F at R1 5' end, 806R at R2 5' end)",
        "fw_primer": PRIMER_515F,
        "rv_primer": PRIMER_806R,
        "r1_search_range": 30,     # Search at 5' end only
        "r2_trunc_len": 150,
        "extract_v4_only": False,  # Already V4 only
    },
    "v45": {
        "description": "V4-V5 amplicons (extract V4 from 515F-926R reads)",
        "fw_primer": PRIMER_515F,
        "rv_primer": PRIMER_926R,
        "r1_search_range": 30,     # Search at 5' end only
        "r2_trunc_len": 150,       # After V4 extraction, same as other amplicons
        "extract_v4_only": True,   # Need to remove V5 region
        "v4_end_primer": PRIMER_806R,      # To find V4 end in R1
        "v4_start_primer": PRIMER_806R_RC, # To find V4 start in R2 (806R rev comp)
    },
}

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
        description="Extract V4 region from paired-end FASTQ by primer search. Supports v34, v4, v45 amplicons."
    )
    parser.add_argument("sample",          help="Sample name (input: sample_R1.fastq.gz / sample_R2.fastq.gz)")
    parser.add_argument("outdir",          help="Output directory")
    parser.add_argument("--amplicon",      choices=["v34", "v4", "v45"], default="v34",
                        help="Amplicon type: v34 (V3-V4), v4 (V4 only), v45 (V4-V5) (default: v34)")
    parser.add_argument("--fw",            default=None,
                        help="Override forward primer (default: auto-selected based on amplicon)")
    parser.add_argument("--rv",            default=None,
                        help="Override reverse primer (default: auto-selected based on amplicon)")
    parser.add_argument("--mismatches",    type=int, default=2,
                        help="Max mismatches allowed per primer (default: 2)")
    parser.add_argument("--r2-rv-window",  type=int, default=30,
                        help="Search window for reverse primer on R2 5' end (default: 30 bp)")
    parser.add_argument("--r2-trunc-len",  type=int, default=None,
                        help="Truncate R2 to this length after primer removal (default: auto-selected based on amplicon)")

    args = parser.parse_args()

    # ── get amplicon configuration ──
    config = AMPLICON_CONFIG[args.amplicon]

    # Use config values, allow CLI overrides
    fw = args.fw if args.fw else config["fw_primer"]
    rv = args.rv if args.rv else config["rv_primer"]
    r1_search_range = config["r1_search_range"]
    r2_trunc = args.r2_trunc_len if args.r2_trunc_len else config["r2_trunc_len"]

    # ── derive paths from sample name ──
    r1_path = f"{args.sample}_R1.fastq.gz"
    r2_path = f"{args.sample}_R2.fastq.gz"

    for p in (r1_path, r2_path):
        if not os.path.isfile(p):
            log.error(f"File not found: {p}")
            sys.exit(1)
    os.makedirs(args.outdir, exist_ok=True)

    log.info("=" * 60)
    log.info(f"Amplicon type     : {args.amplicon.upper()} - {config['description']}")
    log.info("=" * 60)
    log.info(f"Sample            : {args.sample}")
    log.info(f"Forward primer    : {fw}  ({len(fw)} bp)")
    log.info(f"Reverse primer    : {rv}  ({len(rv)} bp)")
    log.info(f"Max mismatches    : {args.mismatches}")
    log.info(f"R1 search range   : {'full read' if r1_search_range is None else f'{r1_search_range} bp (5-prime end)'}")
    log.info(f"R2 RV window      : {args.r2_rv_window} bp")
    log.info(f"R2 truncate len   : {r2_trunc} bp")
    log.info(f"Input             : {r1_path}  /  {r2_path}")

    out_r1 = os.path.join(args.outdir, f"{args.sample}_trimmed_R1.fastq.gz")
    out_r2 = os.path.join(args.outdir, f"{args.sample}_trimmed_R2.fastq.gz")

    # ── V4-V5 specific: need additional primers to extract V4 only ──
    extract_v4_only = config.get("extract_v4_only", False)
    if extract_v4_only:
        v4_end_primer = config["v4_end_primer"]      # 806R - to find V4 end in R1
        v4_start_primer = config["v4_start_primer"]  # 806R_rc - to find V4 start in R2
        log.info(f"V4-V5 mode: will extract V4 region only (remove V5)")
        log.info(f"  R1 V4-end primer  : {v4_end_primer}")
        log.info(f"  R2 V4-start primer: {v4_start_primer}")

    # ── counters ──
    total = kept = 0
    r1_no_fw = r2_no_rv = r1_empty = r2_empty = 0
    r1_no_v4end = r2_no_v4start = 0  # V4-V5 specific counters
    t0 = time.time()

    with open_fq(r1_path) as fh1, open_fq(r2_path) as fh2, \
         gzip.open(out_r1, "wt") as oh1, gzip.open(out_r2, "wt") as oh2:

        for (h1, s1, q1), (h2, s2, q2) in zip(iter_fastq(fh1), iter_fastq(fh2)):
            total += 1

            # ── R1: Search for forward primer, trim if found, use as-is if not ──
            pos_fw = find_primer(s1, fw, args.mismatches, search_range=r1_search_range)
            if pos_fw < 0:
                # Primer not found - assume already trimmed, use R1 as-is
                r1_no_fw += 1
                s1_out = s1
                q1_out = q1
            else:
                trim_start_r1 = pos_fw + len(fw)
                s1_out = s1[trim_start_r1:]
                q1_out = q1[trim_start_r1:]

            # ── V4-V5: Find 806R in R1 and truncate before it (remove V5) ──
            if extract_v4_only and s1_out:
                # Search for 806R starting after a minimum V4 length (~200bp)
                pos_v4end = find_primer(s1_out, v4_end_primer, args.mismatches, search_range=None)
                if pos_v4end > 0:
                    # Truncate before the 806R primer (keep only V4)
                    s1_out = s1_out[:pos_v4end]
                    q1_out = q1_out[:pos_v4end]
                else:
                    # 806R not found - R1 may be too short to contain V5, use as-is
                    r1_no_v4end += 1

            if not s1_out:
                r1_empty += 1
                continue

            # ── R2: Search for reverse primer, trim if found, use as-is if not ──
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

            # ── V4-V5: Find 806R_revcomp in R2 and keep after it (skip V5) ──
            if extract_v4_only and s2_trimmed:
                # Search for 806R_revcomp in the first ~150bp (V5 is ~100bp)
                pos_v4start = find_primer(s2_trimmed, v4_start_primer, args.mismatches, search_range=150)
                if pos_v4start >= 0:
                    # Skip past the 806R_revcomp primer to get V4 only
                    v4_start_pos = pos_v4start + len(v4_start_primer)
                    s2_trimmed = s2_trimmed[v4_start_pos:]
                    q2_trimmed = q2_trimmed[v4_start_pos:]
                else:
                    # 806R_revcomp not found - may be too short, skip this read pair
                    r2_no_v4start += 1
                    continue

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

    log.info("-" * 60)
    log.info(f"  Total pairs processed       : {total:>10,}")
    log.info(f"  R1 - FW primer not found    : {r1_no_fw:>10,}  ({r1_no_fw/total*100:5.2f}%)  (used as-is)")
    log.info(f"  R2 - RV primer not found    : {r2_no_rv:>10,}  ({r2_no_rv/total*100:5.2f}%)  (used as-is)")
    if extract_v4_only:
        log.info(f"  R1 - 806R not found (V4end) : {r1_no_v4end:>10,}  ({r1_no_v4end/total*100:5.2f}%)  (used as-is)")
        log.info(f"  R2 - 806Rrc not found (V4)  : {r2_no_v4start:>10,}  ({r2_no_v4start/total*100:5.2f}%)  (discarded)")
    log.info(f"  R1 - empty after trim       : {r1_empty:>10,}  ({r1_empty/total*100:5.2f}%)")
    log.info(f"  R2 - empty after trim       : {r2_empty:>10,}  ({r2_empty/total*100:5.2f}%)")
    log.info(f"  R2 - truncated to           : {r2_trunc:>10,} bp")
    log.info(f"  Pairs kept                  : {kept:>10,}  ({kept/total*100:5.2f}%)")
    log.info(f"  Output                      : {out_r1}")
    log.info(f"                                {out_r2}")
    log.info(f"  Elapsed                     : {elapsed:>9.1f} s")
    log.info("-" * 60)


if __name__ == "__main__":
    main()
