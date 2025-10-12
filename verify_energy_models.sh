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
COMPILER=${CC:-g++}

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

echo "Comparing develop turner2004 output with main..."
if diff -q "$DEV_T2004_OUT" "$MAIN_OUT" >/dev/null; then
	echo "✔ turner2004 output matches main branch."
else
	echo "✖ turner2004 output differs from main branch." >&2
	diff -u "$MAIN_OUT" "$DEV_T2004_OUT" || true
	exit 1
fi

echo "Comparing develop turner1999 output with legacy..."
if diff -q "$DEV_T1999_OUT" "$LEGACY_OUT" >/dev/null; then
	echo "✔ turner1999 output matches legacy branch."
else
	echo "✖ turner1999 output differs from legacy branch." >&2
	diff -u "$LEGACY_OUT" "$DEV_T1999_OUT" || true
	exit 1
fi

echo "All checks passed. Outputs match both reference branches."
