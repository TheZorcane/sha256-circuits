#!/bin/sh
# Negative controls quoted in the paper (run from the artifact root; reproduce.sh runs this in quick mode).
#  1. twin-cut certificate (96,273 line): 40 mutated lemmas (an operand complement of a twinned gate flipped, or the
#     lemma polarity flipped; twin_cut_certificate/negctl.py, seed 1) must all be refuted by the SAT miter;
#  2. GF(2) checkers: 21 functionally different mutants of the 96,273 circuit (rewrite_checks/mutants.py, seed 1)
#     must be rejected by both lincheck.py and independent/check2.py.
# Prints "NEGATIVE CONTROLS OK" when both hold.  Needs CaDiCaL on PATH; about 1-2 minutes.
set -u
HERE=$(cd "$(dirname "$0")/.." && pwd)
R="$HERE/circuits/proofbase/sha256_compress_96273_proofbase.xag"
N="$HERE/circuits/sha256_compress_96273.xag"
export PYTHONPATH="$HERE/verifier" CERT_REF="$R" CERT_NEW="$N"
cd "$HERE/twin_cut_certificate" || exit 1
rm -rf work; mkdir work
python3 prep.py > /dev/null 2>&1 || { echo "prep.py failed"; exit 1; }
python3 negctl.py 40 1 | tee work/negctl.log
rm -rf work/neg
ok1=0; tail -1 work/negctl.log | grep -q "'not_equivalent': 40, 'other': 0" && ok1=1
rm -rf work
cd "$HERE" || exit 1
python3 rewrite_checks/mutants.py "$R" "$N" 21 1 | tee /dev/stderr | tail -1 | grep -q "MUTANTS REJECTED: 21 of 21" && ok2=1 || ok2=0
echo "twin-cut negative controls: $([ $ok1 = 1 ] && echo '40 of 40 refuted' || echo FAIL)"
echo "GF(2) mutants: $([ $ok2 = 1 ] && echo '21 of 21 rejected by both checkers' || echo FAIL)"
[ $ok1 = 1 ] && [ $ok2 = 1 ] && echo "NEGATIVE CONTROLS OK"
