#!/bin/bash
# Re-run every SAW proof that links the verifier's SHA-256 reference to Galois's Cryptol FIPS 180-4 spec.
#   usage: CRYPTOL_SPECS=/path/to/cryptol-specs SAW=/path/to/saw ./run.sh
# Needs: SAW 1.6 (with Cryptol 3.6), yices 2.8 on PATH, python3, and a checkout of
# https://github.com/GaloisInc/cryptol-specs at commit 8b74f5fc2b77fc7bbb587ba2217990b1963e9d6c.
# SAW_JOBS (default 3) sets how many SAW processes run at once (each runs one single-threaded yices);
# about 10 minutes wall with 3 jobs; about 25 minutes with SAW_JOBS=1 (sum of the per-run wall times).
# Exit status: 0 only if the cryptol-specs commit is the one above, every SAW run exits 0, the negative
# control rejects the mutant, all 329 reference links hold, and summarize.py finds every required log with
# the expected theorems (see summarize.py).  Anything else gives exit status 1.
set -u -o pipefail
FAILED=0
SPEC_COMMIT=8b74f5fc2b77fc7bbb587ba2217990b1963e9d6c
V="$(cd "$(dirname "$0")" && pwd)"
B="$V/build"; L="${LOGS:-$V/out}"   # build/ must stay inside this folder (compose.saw includes ../saw/...)
: "${CRYPTOL_SPECS:?set CRYPTOL_SPECS to a cryptol-specs checkout}"
export SAW="${SAW:-saw}"
SAW_JOBS="${SAW_JOBS:-3}"
mkdir -p "$B" "$L"
export CRYPTOLPATH="$CRYPTOL_SPECS:$V/cry:$B:$B/plain:$B/tail:$B/split"
# never summarize logs left over from an earlier run
rm -f "$L"/frags_[SI]_*.log "$L"/compose_*.log "$L"/negative.log "$L"/ref_link.log "$L"/gen.log "$L"/env.txt "$L"/summary.txt

COMMIT="$(cd "$CRYPTOL_SPECS" && git rev-parse HEAD 2>/dev/null || echo unknown)"
{ echo "cryptol-specs commit: $COMMIT"
  "$SAW" --version 2>&1
  yices --version | head -1; } | tee "$L/env.txt"
if [ "$COMMIT" != "$SPEC_COMMIT" ]; then
  echo "cryptol-specs checkout is at $COMMIT, expected $SPEC_COMMIT: FAILED"; exit 1
fi

echo "== generate fragments and AIGER files (each validated by simulation), Cryptol lemma modules"
python3 "$V/gen.py" "$B" > "$L/gen.log" || { echo "gen.py FAILED"; exit 1; }
grep -q '"identical_to_inlined_fragments_plain": true' "$L/gen.log" || { echo "inline identity check FAILED"; exit 1; }

for var in plain tail split; do
  printf 'import "Lemmas_%s_S.cry" as S;\nimport "Lemmas_%s_I.cry" as I;\nlet var = "%s";\ninclude "%s/saw/compose.saw";\n' \
    "$var" "$var" "$var" "$V" > "$B/compose_$var.saw"
done
cp "$V/saw/negative.saw" "$B/negative.saw"

run() { python3 "$V/run_saw.py" "$@"; }
PIDS=""
waitall() {  # wait for each background SAW run; any nonzero status is a failure
  for p in $PIDS; do wait "$p" || { echo "a SAW run (pid $p) FAILED"; FAILED=1; }; done
  PIDS=""
}
for var in plain tail split; do
  for st in S I; do
    run 3600 "$L/frags_${st}_$var.log" "$B/$var/frags_$st.saw" & PIDS="$PIDS $!"
    while [ "$(jobs -rp | wc -l)" -ge "$SAW_JOBS" ]; do sleep 1; done
  done
done
waitall
for var in plain tail split; do
  run 3600 "$L/compose_$var.log" "$B/compose_$var.saw" & PIDS="$PIDS $!"
  while [ "$(jobs -rp | wc -l)" -ge "$SAW_JOBS" ]; do sleep 1; done
done
waitall
run 600 "$L/negative.log" "$B/negative.saw" || { echo "negative control FAILED"; FAILED=1; }

echo "== SAT reference <-> SAW fragment link (329 lemmas of the three circuits; AIGER files of $B)"
python3 "$V/check_ref_link.py" all "$B" > "$L/ref_link.log" 2>&1 || { echo "reference link FAILED"; FAILED=1; }
tail -1 "$L/ref_link.log"

echo "== summary (summarize.py checks every required log and theorem count)"
python3 "$V/summarize.py" "$L" > "$L/summary.txt" 2>&1 || FAILED=1
cat "$L/summary.txt"
if [ "$FAILED" -eq 0 ]; then echo "RUN.SH: ALL REQUIRED PROOFS SUCCEEDED"; else echo "RUN.SH: FAILED"; fi
exit "$FAILED"
