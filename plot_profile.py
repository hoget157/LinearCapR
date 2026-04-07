#!/usr/bin/env python3
"""Plot LinearCapR structural-context profiles.

Two common usage modes are supported:

1. Plot from an existing LinearCapR profile file:
   python plot_profile.py --profile example.profile --output-prefix example

2. Run LinearCapR from a FASTA file and then plot:
   python plot_profile.py --fasta example.fa --output-prefix example

An optional reference secondary structure in dot-bracket format can be supplied
to draw a one-line context-label track above the profile plot.
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

matplotlib.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": [
            "Times New Roman",
            "Times New Roman PS MT",
            "Times",
            "Nimbus Roman No9 L",
            "DejaVu Serif",
        ],
        "font.size": 12,
        "axes.titlesize": 15,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "legend.title_fontsize": 12,
        "figure.titlesize": 16,
    }
)


CONTEXT_ORDER = [
    ("Hairpin", "H", "#1b9e77"),
    ("Bulge", "B", "#d95f02"),
    ("Internal", "I", "#7570b3"),
    ("Multiloop", "M", "#e7298a"),
    ("Exterior", "E", "#66a61e"),
    ("Stem", "S", "#1f78b4"),
]

SHORT_TO_COLOR = {short: color for _, short, color in CONTEXT_ORDER}
SHORT_TO_INDEX = {short: idx for idx, (_, short, _) in enumerate(CONTEXT_ORDER)}


class StructureParser:
    def __init__(self, structure: str):
        self.structure = structure
        self.n = len(structure)
        self.pt = self._make_pair_table()

    def _make_pair_table(self) -> list[int]:
        pt = [0] * (self.n + 1)
        stack: list[int] = []
        for i, char in enumerate(self.structure, 1):
            if char == "(":
                stack.append(i)
            elif char == ")" and stack:
                j = stack.pop()
                pt[i] = j
                pt[j] = i
        return pt

    def get_pairs(self) -> list[tuple[int, int]]:
        return [(i, j) for i, j in enumerate(self.pt[1:], 1) if j > i]

    def is_hairpin_closing(self, i: int, j: int) -> bool:
        if self.pt[i] != j or i >= j:
            return False
        return all(self.pt[k] == 0 for k in range(i + 1, j))

    def is_multiloop_start(self, i: int, j: int) -> bool:
        if self.pt[i] != j or i >= j:
            return False
        branches = 0
        k = i + 1
        while k < j:
            if self.pt[k] > k:
                branches += 1
                k = self.pt[k] + 1
            else:
                k += 1
        return branches >= 2

    def is_internal_loop(self, i: int, j: int, h: int, l: int) -> bool:
        if not (self.pt[i] == j and self.pt[h] == l):
            return False
        if not (i < h < l < j):
            return False
        if h == i + 1 and l == j - 1:
            return True
        if any(self.pt[k] != 0 for k in range(i + 1, h)):
            return False
        if any(self.pt[k] != 0 for k in range(l + 1, j)):
            return False
        return True


class StructureContextLabeler(StructureParser):
    def __init__(self, structure: str):
        super().__init__(structure)
        self.pairs = self.get_pairs()
        self.hairpin_pairs: list[tuple[int, int]] = []
        self.multiloop_starts: list[tuple[int, int]] = []
        self.internal_loops: list[tuple[int, int, int, int]] = []
        self.outermost_pairs: list[tuple[int, int]] | None = None
        self._analyze()

    def _analyze(self) -> None:
        for i, j in self.pairs:
            if self.is_hairpin_closing(i, j):
                self.hairpin_pairs.append((i, j))
            elif self.is_multiloop_start(i, j):
                self.multiloop_starts.append((i, j))
            for h, l in self.pairs:
                if h > i and l < j and self.is_internal_loop(i, j, h, l) and not (h == i + 1 and l == j - 1):
                    self.internal_loops.append((i, j, h, l))

    def is_stem(self, k: int) -> bool:
        return self.pt[k] != 0

    def is_in_hairpin(self, k: int) -> bool:
        return any(i < k < j for i, j in self.hairpin_pairs)

    def is_in_bulge(self, k: int) -> bool:
        for i, j, h, l in self.internal_loops:
            left_unpaired = list(range(i + 1, h))
            right_unpaired = list(range(l + 1, j))
            is_bulge = (len(left_unpaired) > 0) != (len(right_unpaired) > 0)
            if is_bulge and (k in left_unpaired or k in right_unpaired):
                return True
        return False

    def is_in_internal(self, k: int) -> bool:
        for i, j, h, l in self.internal_loops:
            left_unpaired = list(range(i + 1, h))
            right_unpaired = list(range(l + 1, j))
            if left_unpaired and right_unpaired and (k in left_unpaired or k in right_unpaired):
                return True
        return False

    def is_in_exterior(self, k: int) -> bool:
        if self.outermost_pairs is None:
            self.outermost_pairs = []
            for i, j in self.pairs:
                is_outermost = True
                for ii, jj in self.pairs:
                    if ii < i and jj > j:
                        is_outermost = False
                        break
                if is_outermost:
                    self.outermost_pairs.append((i, j))
        return not any(i <= k <= j for i, j in self.outermost_pairs)

    def is_in_multiloop(self, k: int) -> bool:
        for i, j in self.multiloop_starts:
            if not (i < k < j):
                continue
            branch_pairs = []
            for h, l in self.pairs:
                if i < h and l < j:
                    is_branch = True
                    for ii, jj in self.pairs:
                        if i < ii < h and l < jj < j:
                            is_branch = False
                            break
                    if is_branch:
                        branch_pairs.append((h, l))
            if not any(h <= k <= l for h, l in branch_pairs):
                return True
        return False

    def get_all_labels(self) -> list[str]:
        labels = []
        for k in range(1, self.n + 1):
            if self.is_stem(k):
                labels.append("S")
            elif self.is_in_hairpin(k):
                labels.append("H")
            elif self.is_in_bulge(k):
                labels.append("B")
            elif self.is_in_internal(k):
                labels.append("I")
            elif self.is_in_multiloop(k):
                labels.append("M")
            elif self.is_in_exterior(k):
                labels.append("E")
            else:
                labels.append("E")
        return labels


def parse_profile(path: Path) -> tuple[str, dict[str, list[float]]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Unexpected profile format: {path}")

    sequence_name = lines[0][1:]
    values: dict[str, list[float]] = {}
    label_map = {"Multibranch": "Multiloop"}
    for row in lines[1:]:
        parts = row.split()
        if len(parts) < 2:
            continue
        context = label_map.get(parts[0], parts[0])
        values[context] = [float(x) for x in parts[1:]]
    return sequence_name, values


def read_dot_bracket(path: Path) -> str:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) == 1:
        return lines[0]
    if len(lines) >= 3 and lines[0].startswith(">"):
        return lines[2]
    return lines[-1]


def build_xticks(length: int, step: int = 20) -> list[int]:
    ticks = [1]
    cur = step
    while cur < length:
        ticks.append(cur)
        cur += step
    if ticks[-1] != length:
        ticks.append(length)
    return ticks


def run_lincapr(
    fasta: Path,
    output_profile: Path,
    beam_size: int,
    energy_model: str,
    lincapr_bin: Path,
) -> None:
    cmd = [str(lincapr_bin), str(fasta), str(output_profile), str(beam_size), "--energy", energy_model]
    subprocess.run(cmd, check=True)


def plot_overlay(
    sequence_name: str,
    values: dict[str, list[float]],
    output_prefix: Path,
    figure_title: str,
    reference_track_name: str | None,
    reference_labels: list[str] | None,
) -> None:
    lengths = {len(v) for v in values.values()}
    if len(lengths) != 1:
        raise ValueError("Context rows have different lengths")
    n = lengths.pop()
    x = list(range(1, n + 1))

    n_tracks = int(reference_labels is not None)
    height_ratios = [0.32] * n_tracks + [3.2]
    fig, axes = plt.subplots(
        n_tracks + 1,
        1,
        figsize=(13, 6.8 if n_tracks else 5.8),
        sharex=True,
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.05},
    )
    if hasattr(axes, "ravel"):
        axes = list(axes.ravel())
    else:
        axes = [axes]

    cmap = ListedColormap([color for _, _, color in CONTEXT_ORDER])
    track_axes = axes[:-1]
    main_ax = axes[-1]

    if reference_labels is not None:
        idx_row = [[SHORT_TO_INDEX[ch] for ch in reference_labels]]
        ax = track_axes[0]
        ax.imshow(
            idx_row,
            aspect="auto",
            cmap=cmap,
            interpolation="nearest",
            extent=(0.5, n + 0.5, 0, 1),
        )
        ax.set_yticks([])
        ax.set_ylabel(reference_track_name or "Reference labels", rotation=0, labelpad=80, va="center", fontsize=11)
        for spine in ax.spines.values():
            spine.set_visible(False)

    for context_name, short_name, color in CONTEXT_ORDER:
        y = values.get(context_name)
        if y is None:
            raise ValueError(f"Missing context row: {context_name}")
        main_ax.plot(x, y, color=color, linewidth=1.8, label=short_name)

    main_ax.set_ylim(0.0, 1.0)
    main_ax.set_xlim(1, n)
    main_ax.set_ylabel("Probability", fontsize=13)
    main_ax.set_xlabel("Position (1-origin)", fontsize=13)
    main_ax.set_xticks(build_xticks(n))
    main_ax.grid(True, alpha=0.25, linewidth=0.6)

    patch_handles = [mpatches.Patch(color=color, label=short) for _, short, color in CONTEXT_ORDER]
    main_ax.legend(handles=patch_handles, title="Contexts", loc="upper right", ncol=3, frameon=True)

    fig.suptitle(figure_title or f"LinearCapR profile: {sequence_name}")
    fig.tight_layout(rect=(0.08, 0, 1, 0.965))
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_prefix) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(output_prefix) + ".pdf", bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, help="Existing LinearCapR profile file")
    parser.add_argument("--fasta", type=Path, help="FASTA input to run through LinearCapR before plotting")
    parser.add_argument("--output-prefix", type=Path, required=True, help="Output path prefix (without extension)")
    parser.add_argument("--beam-size", type=int, default=100, help="Beam size for --fasta mode (default: 100)")
    parser.add_argument(
        "--energy",
        choices=["turner2004", "turner1999"],
        default="turner2004",
        help="Energy model for --fasta mode (default: turner2004)",
    )
    parser.add_argument(
        "--lincapr-bin",
        type=Path,
        default=Path("./LinCapR"),
        help="Path to the LinearCapR executable for --fasta mode",
    )
    parser.add_argument(
        "--reference-secstruct",
        type=Path,
        default=None,
        help="Optional reference secondary structure in dot-bracket format",
    )
    parser.add_argument(
        "--reference-track-name",
        default="Reference structure labels",
        help="Display name for the optional reference label track",
    )
    parser.add_argument(
        "--figure-title",
        default="LinearCapR profile",
        help="Figure title",
    )
    args = parser.parse_args()
    if args.profile is None and args.fasta is None:
        parser.error("one of --profile or --fasta is required")
    return args


def main() -> None:
    args = parse_args()

    if args.profile is not None:
        profile_path = args.profile
    else:
        args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="lincapr-plot-") as tmpdir:
            tmp_profile = Path(tmpdir) / "tmp.profile"
            run_lincapr(args.fasta, tmp_profile, args.beam_size, args.energy, args.lincapr_bin)
            sequence_name, values = parse_profile(tmp_profile)
            reference_labels = None
            if args.reference_secstruct is not None:
                structure = read_dot_bracket(args.reference_secstruct)
                reference_labels = StructureContextLabeler(structure).get_all_labels()
                if len(reference_labels) != len(next(iter(values.values()))):
                    raise ValueError("Reference structure length does not match the profile length")
            plot_overlay(
                sequence_name,
                values,
                args.output_prefix,
                args.figure_title,
                args.reference_track_name,
                reference_labels,
            )
            return

    sequence_name, values = parse_profile(profile_path)
    reference_labels = None
    if args.reference_secstruct is not None:
        structure = read_dot_bracket(args.reference_secstruct)
        reference_labels = StructureContextLabeler(structure).get_all_labels()
        if len(reference_labels) != len(next(iter(values.values()))):
            raise ValueError("Reference structure length does not match the profile length")

    plot_overlay(
        sequence_name,
        values,
        args.output_prefix,
        args.figure_title,
        args.reference_track_name,
        reference_labels,
    )


if __name__ == "__main__":
    main()
