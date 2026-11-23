#!/bin/sh
# usage: run_certificate.sh REF.xag NEW.xag OUTDIR [JOBS]
# Twin-cut LRAT certificate that NEW computes the same function as REF, where NEW keeps every AND function of
# REF and differs only in its linear (XOR) layer.  One SAT lemma per twin node (verifier/equiv.py: miter,
# CaDiCaL LRAT, strict RUP check + second checker), then compose.py re-derives the induction with its own
# circuit builder.  Keeps logs + JSON manifests in OUTDIR, deletes CNF/LRAT after composition.
# Needs about 1 GB of free disk while running; 3-6 minutes on a laptop with 4 jobs.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); N=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
mkdir -p "$3"; O=$(cd "$3" && pwd); J=${4:-4}
export CERT_REF="$R" CERT_NEW="$N" PYTHONPATH="$HERE/../verifier"
cd "$HERE"; rm -rf work; mkdir work
if command -v sha256sum >/dev/null 2>&1; then sha256sum "$R" "$N"; else shasum -a 256 "$R" "$N"; fi > "$O/input_sha256.txt"
python3 prep.py > "$O/prep.log" 2>&1
python3 run.py work/out --jobs "$J" > "$O/run.log" 2>&1
if python3 compose.py work/out > "$O/compose.log" 2>&1; then echo compose_exit=0 >> "$O/compose.log"; else echo compose_exit=FAIL >> "$O/compose.log"; exit 1; fi
cp work/plan.json "$O/"
(cd work/out && find . -name '*.json' | tar czf "$O/manifests.tgz" -T -)
rm -rf work
tail -2 "$O/run.log"; tail -6 "$O/compose.log"
