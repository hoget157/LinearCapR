#!/usr/bin/env python3
"""
Compare LinearCapR output profiles.

Usage:
    ./compare_profiles.py file_a file_b

The script parses the profile text files produced by LinearCapR (lines beginning
with profile labels such as "Bulge", "Hairpin", etc.), aligns sequences by their
headers, and reports statistics for each profile type:
    * maximum absolute difference
    * mean absolute difference
    * RMS difference
    * indices where the maximum occurs

Any structural profile not present in both files is flagged. The script exits
with status 0 if all sequences and profile entries match bit‑exactly, otherwise
it returns 1.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

ProfileMap = Dict[str, Dict[str, List[float]]]


def parse_profile(path: Path) -> ProfileMap:
    """Parse a LinearCapR output file into a nested dict."""
    current_seq: str | None = None
    profiles: ProfileMap = {}

    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_seq = line[1:]
                profiles[current_seq] = {}
                continue
            if current_seq is None:
                raise ValueError(f"Encountered data before sequence header in {path}: {line}")
            label, *values = line.split()
            try:
                profiles[current_seq][label] = [float(v) for v in values]
            except ValueError as exc:
                raise ValueError(f"Failed to parse floats for {label} in {path}") from exc
    return profiles


def compare_vectors(vec_a: List[float], vec_b: List[float]) -> Tuple[float, float, float, List[int]]:
    if len(vec_a) != len(vec_b):
        raise ValueError(f"Length mismatch: {len(vec_a)} vs {len(vec_b)}")

    deltas = [abs(a - b) for a, b in zip(vec_a, vec_b)]
    if not deltas:
        return 0.0, 0.0, 0.0, []

    max_delta = max(deltas)
    max_indices = [i for i, d in enumerate(deltas) if d == max_delta]
    mean_delta = sum(deltas) / len(deltas)
    rms_delta = math.sqrt(sum(d * d for d in deltas) / len(deltas))
    return max_delta, mean_delta, rms_delta, max_indices


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare LinearCapR profile outputs.")
    parser.add_argument("file_a", type=Path)
    parser.add_argument("file_b", type=Path)
    parser.add_argument("--tolerance", "-t", type=float, default=1e-20, help="Allowed max absolute deviation (default: 1e-5).")
    parser.add_argument("--verbose", "-v", action="store_true", help="Show per-index differences when mismatches exceed tolerance.")
    args = parser.parse_args()

    profiles_a = parse_profile(args.file_a)
    profiles_b = parse_profile(args.file_b)

    exit_code = 0

    keys_a = set(profiles_a)
    keys_b = set(profiles_b)

    only_a = sorted(keys_a - keys_b)
    only_b = sorted(keys_b - keys_a)

    if only_a:
        print(f"Sequences only in {args.file_a}: {', '.join(only_a)}")
        exit_code = 1
    if only_b:
        print(f"Sequences only in {args.file_b}: {', '.join(only_b)}")
        exit_code = 1

    for seq in sorted(keys_a & keys_b):
        prof_a = profiles_a[seq]
        prof_b = profiles_b[seq]
        labels = set(prof_a) | set(prof_b)

        for label in sorted(labels):
            if label not in prof_a:
                print(f"[{seq}] Missing label {label} in {args.file_a}")
                exit_code = 1
                continue
            if label not in prof_b:
                print(f"[{seq}] Missing label {label} in {args.file_b}")
                exit_code = 1
                continue

            vec_a = prof_a[label]
            vec_b = prof_b[label]

            if len(vec_a) != len(vec_b):
                print(f"[{seq}][{label}] length mismatch: {len(vec_a)} vs {len(vec_b)}")
                exit_code = 1
                continue

            max_delta, mean_delta, rms_delta, max_indices = compare_vectors(vec_a, vec_b)
            if max_delta <= args.tolerance:
                continue

            exit_code = 1
            print(
                f"[{seq}][{label}] max Δ={max_delta:.6g} "
                f"(mean={mean_delta:.6g}, rms={rms_delta:.6g}) exceeds tolerance {args.tolerance:g} "
                f"at indices {max_indices}"
            )
            if args.verbose:
                for idx in max_indices:
                    print(f"    idx {idx}: {vec_a[idx]} vs {vec_b[idx]}")

    if exit_code == 0:
        print("All sequences and profile values match exactly.")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
