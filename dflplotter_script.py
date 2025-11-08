#!/usr/bin/env python3
"""
Plot DFLpred scores along a protein sequence, call low-flexibility (DFL) regions,
and optionally compare "case-annotated" DFL calls (lowercase = DFL).

Inputs:
  - protein ID (string)
  - sequence (string; lowercase letters indicate DFL by annotation)
  - scores (comma-separated string or a file path with one score per residue)
Outputs:
  - PNG plot
  - Text summary of called regions (stdout)

"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt

THRESHOLD_DEFAULT = 0.18

def parse_scores(scores_arg: str):
    p = Path(scores_arg)
    if p.exists():
        with p.open() as fh:
            vals = [line.strip() for line in fh]
        flat = ",".join(vals)
    else:
        flat = scores_arg
    scores = [float(x.strip()) for x in flat.split(",") if x.strip() != ""]
    return scores

def call_regions_below_threshold(flags):
    regions = []
    start = None
    n = len(flags)
    for i, flag in enumerate(flags, start=1):  # 1-based indexing
        if flag and start is None:
            start = i
        if (not flag or i == n) and start is not None:
            end = i if flag and i == n else i - 1
            regions.append((start, end))
            start = None
    return regions

def main():
    ap = argparse.ArgumentParser(
        description="Plot DFLpred scores and call DFL regions (scores < threshold)."
    )
    ap.add_argument("--protein-id", required=True, help="Protein identifier (string).")
    ap.add_argument("--sequence", required=True,
                    help="AAs (letters). Lowercase indicates DFL by annotation.")
    ap.add_argument("--scores", required=True,
                    help="Comma-separated scores or a file path with one score per line.")
    ap.add_argument("--threshold", type=float, default=THRESHOLD_DEFAULT,
                    help=f"DFL score threshold (default: {THRESHOLD_DEFAULT}).")
    ap.add_argument("--out", default="dflpred_plot.png",
                    help="Output plot filename (PNG). Default: dflpred_plot.png")
    args = ap.parse_args()

    protein_id = args.protein_id
    seq = args.sequence.strip()
    scores = parse_scores(args.scores)

    if len(seq) != len(scores):
        raise SystemExit(
            f"Error: sequence length ({len(seq)}) != number of scores ({len(scores)})."
        )

    # Calls from annotation (lowercase) and threshold
    is_dfl_case = [aa.islower() for aa in seq]
    is_dfl_thr = [s < args.threshold for s in scores]

    # Consistency check
    agree = sum(int(a == b) for a, b in zip(is_dfl_case, is_dfl_thr))
    if agree != len(seq):
        print(f"Note: {len(seq) - agree} residues differ between case-annotation and threshold calls.")

    # Build contiguous regions using threshold calls
    regions = call_regions_below_threshold(is_dfl_thr)

    # Print summary
    print(protein_id)
    print(f"Length: {len(seq)} aa")
    print(f"DFL threshold: {args.threshold}")
    print("Predicted DFL regions (threshold-based):")
    if regions:
        for (a, b) in regions:
            frag = seq[a-1:b]
            print(f"  {a:>4}-{b:<4}  len={b-a+1:<3}  seq={frag}")
    else:
        print("  none")

    # Plot
    x = list(range(1, len(scores) + 1))
    plt.figure(figsize=(10, 4.5))
    plt.plot(x, scores, linewidth=1.5)
    plt.axhline(args.threshold, linestyle="--", linewidth=1)

    # Shade predicted DFL regions (below threshold)
    for a, b in regions:
        plt.axvspan(a, b, alpha=0.2)

    title_id = protein_id.lstrip(">")
    plt.title(f"DFLpred scores — {title_id}")
    plt.xlabel("Residue number")
    plt.ylabel("DFLpred score")
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    print(f"Saved figure: {args.out}")

if __name__ == "__main__":
    main()


