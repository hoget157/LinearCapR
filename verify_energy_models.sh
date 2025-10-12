#!/usr/bin/env bash
set -euo pipefail

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
	echo "This script must be run inside the LinearCapR git repository." >&2
	exit 1
fi

ROOT="$(git rev-parse --show-toplevel)"
if [[ -n "$(git status --porcelain)" ]]; then
	echo "Error: please start from a clean worktree." >&2
	exit 1
fi

INPUT_DEFAULT="test.fa"
INPUT_REL="${1:-$INPUT_DEFAULT}"
BEAM_SIZE="${2:-100}"

INPUT_PATH="$(python - "$INPUT_REL" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
)"

if [[ ! -f "$INPUT_PATH" ]]; then
	echo "Error: input file '$INPUT_PATH' not found." >&2
	exit 1
fi

TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/lincapr-verify.XXXXXXXX")"
cleanup(){
	if [[ -d "$TMP_ROOT" ]]; then
		if git -C "$ROOT" worktree list --porcelain | grep -q "worktree $TMP_ROOT/main"; then
			git -C "$ROOT" worktree remove --force "$TMP_ROOT/main"
		fi
		if git -C "$ROOT" worktree list --porcelain | grep -q "worktree $TMP_ROOT/legacy"; then
			git -C "$ROOT" worktree remove --force "$TMP_ROOT/legacy"
		fi
		rm -rf "$TMP_ROOT"
	fi
}
trap cleanup EXIT

git -C "$ROOT" worktree add --quiet "$TMP_ROOT/main" main
git -C "$ROOT" worktree add --quiet "$TMP_ROOT/legacy" legacy

MAKE_CMD=${MAKE:-make}
COMPILER=${CC:-clang++}

build_branch(){
	local dir="$1"
	local label="$2"
	echo "Building ${label} branch..."
	("$MAKE_CMD" -C "$dir" clean CC="$COMPILER" >/dev/null 2>&1 || true)
	"$MAKE_CMD" -C "$dir" CC="$COMPILER"
}

build_branch "$ROOT" "develop (current)"
build_branch "$TMP_ROOT/main" "main"
build_branch "$TMP_ROOT/legacy" "legacy"

DEV_T2004_OUT="$TMP_ROOT/develop_turner2004.out"
DEV_T1999_OUT="$TMP_ROOT/develop_turner1999.out"
MAIN_OUT="$TMP_ROOT/main_turner2004.out"
LEGACY_OUT="$TMP_ROOT/legacy_turner1999.out"

echo "Running develop with turner2004 parameters..."
"$ROOT/LinCapR" "$INPUT_PATH" "$DEV_T2004_OUT" "$BEAM_SIZE" --energy turner2004

echo "Running develop with turner1999 parameters..."
"$ROOT/LinCapR" "$INPUT_PATH" "$DEV_T1999_OUT" "$BEAM_SIZE" --energy turner1999

echo "Running main reference..."
"$TMP_ROOT/main/LinCapR" "$INPUT_PATH" "$MAIN_OUT" "$BEAM_SIZE"

echo "Running legacy reference..."
"$TMP_ROOT/legacy/LinearCapR" "$INPUT_PATH" "$LEGACY_OUT" "$BEAM_SIZE"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPARE_PY="$SCRIPT_DIR/compare_profiles.py"

overall_status=0

run_compare(){
	local description="$1"
	local file_a="$2"
	local file_b="$3"

	echo "Comparing $description..."
	if [[ -x "$COMPARE_PY" ]]; then
		if "$COMPARE_PY" "$file_a" "$file_b"; then
			echo "✔ $description within tolerance."
		else
			echo "✖ $description differs beyond tolerance." >&2
			overall_status=1
		fi
	else
		echo "Warning: compare_profiles.py not found or not executable; falling back to diff." >&2
		if diff -q "$file_a" "$file_b" >/dev/null; then
			echo "✔ $description matches exactly."
		else
			echo "✖ $description differs (exact diff not shown)." >&2
			overall_status=1
		fi
	fi
	echo
}

run_compare "develop turner2004 vs main" "$DEV_T2004_OUT" "$MAIN_OUT"
run_compare "develop turner1999 vs legacy" "$DEV_T1999_OUT" "$LEGACY_OUT"

if [[ $overall_status -eq 0 ]]; then
	echo "All checks passed. Outputs match within tolerance."
else
	echo "Differences detected. See reports above." >&2
	exit 1
fi
