#!/bin/sh
# usage: walk.sh BASE.xag PROVEN.xag OUT.xag [iterations_per_round rounds]
# Residual linear-layer walk: keeps every AND gate (same AND functions, possibly re-chosen operand pairs inside
# each AND's plane) and re-synthesises the XOR layer.  Census of weight-3/4 linear dependencies -> augmented
# circuit + AND planes -> deterministic Steiner-style local search -> emit -> canonicalize -> two independent
# GF(2) checks against PROVEN (the circuit whose equivalence to SHA-256 is proved).  Deterministic (fixed seeds).
# Default 8 rounds x 5e7 moves: about 5-8 minutes and 2 threads on a laptop.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
SUP="$HERE/.."
export PYTHONPATH="$SUP/verifier"
abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
B=$(abs "$1"); PB=$(abs "$2"); OUT=$(abs "$3"); IT=${4:-50000000}; RO=${5:-8}
BIN="$HERE/bin"; mkdir -p "$BIN"
[ -x "$BIN/deps" ] || cc -O2 -o "$BIN/deps" "$HERE/deps.c"
[ -x "$BIN/steiner" ] || cc -O2 -o "$BIN/steiner" "$HERE/steiner.c" -lm
S=$(mktemp -d "${TMPDIR:-/tmp}/walk.XXXXXX")
N=$(grep -c '^[AX] ' "$B")
"$BIN/deps" "$B" 1 0 > "$S/d.tri" 2> "$S/d.d3log"
for k in 0 1; do "$BIN/deps" "$B" 1 $((64*k+64)) $((64*k)) skip 2> "$S/q4log$k" | grep '^Q' > "$S/q4part$k" & done; wait
for k in 2 3; do "$BIN/deps" "$B" 1 $((64*k+64)) $((64*k)) skip 2> "$S/q4log$k" | grep '^Q' > "$S/q4part$k" & done; wait
cat "$S"/q4part? > "$S/d.q4"
python3 "$HERE/augment.py" "$B" "$S/d.q4" "$S/u.xag" "$S/u.planes"
"$BIN/deps" "$S/u.xag" 1 0 > "$S/u.tri" 2>/dev/null
PLANES="$S/u.planes" "$BIN/steiner" "$S/u.xag" "$N" "$S/u.tri" "$S/u.state" 1 "$IT" "$RO" 2>&1 | tail -1
python3 "$HERE/emit.py" "$S/u.xag" "$S/u.state" "$S/a.xag"
python3 "$HERE/canon.py" "$S/a.xag" "$OUT"
python3 "$SUP/rewrite_checks/lincheck.py" "$PB" "$OUT" | tail -1
(cd "$SUP/rewrite_checks/independent" && python3 check2.py "$PB" "$OUT" | tail -1)
rm -rf "$S"
