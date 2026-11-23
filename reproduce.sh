#!/bin/sh
# Reproduce the results of "Smaller verified Boolean circuits for SHA-256".
#   sh reproduce.sh           quick mode: about 7 minutes on a laptop (up to 4 threads, < 2 GB disk)
#   sh reproduce.sh full      also: rebuild every circuit from the generators and compare bytes (~30 min),
#                             and run the complete compositional proofs including the two fused-round tail
#                             lemmas (about 1 hour and 20 GB of free disk per proof base; 100,658: ~15 min, 6 GB)
#   PUBLIC_RAW=DIR sh reproduce.sh   additionally score the public comparator files in DIR
#                             (download them with comparators/fetch_public_baselines.sh DIR)
# Requirements: Python 3 (standard library only), a C compiler, CaDiCaL 3.x on PATH as `cadical`.
# JOBS (default 4) sets the number of parallel SAT jobs.
set -u
MODE=${1:-quick}
JOBS=${JOBS:-4}
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
OUT="$HERE/out"; mkdir -p "$OUT"
export PYTHONPATH="$HERE/verifier" PYTHONDONTWRITEBYTECODE=1
SUMMARY="$OUT/summary.txt"; : > "$SUMMARY"
step() {  # step NAME EXPECTED_PATTERN COMMAND...
  name=$1; pat=$2; shift 2
  t0=$(date +%s)
  if "$@" > "$OUT/$name.log" 2>&1 && grep -q -- "$pat" "$OUT/$name.log"; then r=PASS; else r=FAIL; fi
  t1=$(date +%s)
  printf '%-34s %s  (%ss)\n' "$name" "$r" $((t1 - t0)) | tee -a "$SUMMARY"
}
sha() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$@"; else shasum -a 256 "$@"; fi; }
shacheck() { if command -v sha256sum >/dev/null 2>&1; then sha256sum -c "$@"; else shasum -a 256 -c "$@"; fi; }

cc -O2 -o verifier/vendor/lrat-check verifier/vendor/lrat-check.c 2>/dev/null
cc -O2 -o verifier/vendor/drat-trim  verifier/vendor/drat-trim.c 2>/dev/null
command -v cadical >/dev/null || echo "WARNING: cadical not on PATH; SAT steps will fail"

# 1. integrity, counts, simulation, Bristol Fashion files
step sha256sums ": OK" shacheck SHA256SUMS
step counts "COUNTS OK" python3 tools/check_circuits.py count
step simulation "SIM OK" python3 tools/check_circuits.py sim 100000
step bristol_reread "BRISTOL OK" python3 tools/check_circuits.py bristol-check

# 2. algebraic checks of the linear-layer rewrite (final circuit vs its proof base)
lin2() {  # four checks (two checkers x two lines); each must end with a line starting with PASS
  n=0
  for c in 96273 96324; do
    r1=$(python3 rewrite_checks/lincheck.py circuits/proofbase/sha256_compress_${c}_proofbase.xag circuits/sha256_compress_${c}.xag | tail -1)
    echo "lincheck $c: $r1"; case "$r1" in PASS*) n=$((n + 1));; esac
    r2=$(cd rewrite_checks/independent && python3 check2.py ../../circuits/proofbase/sha256_compress_${c}_proofbase.xag ../../circuits/sha256_compress_${c}.xag | tail -1)
    echo "check2 $c: $r2"; case "$r2" in PASS*) n=$((n + 1));; esac
  done
  echo "$n of 4 checks passed"; [ "$n" -eq 4 ] && echo "GF2 ALL PASS"
}
step gf2_checks "GF2 ALL PASS" lin2

# 3. twin-cut LRAT certificates: final circuit == proof base
step twin_cut_96273 "compose_exit=0" sh twin_cut_certificate/run_certificate.sh \
     circuits/proofbase/sha256_compress_96273_proofbase.xag circuits/sha256_compress_96273.xag "$OUT/twin_cut_96273" "$JOBS"
step twin_cut_96324 "compose_exit=0" sh twin_cut_certificate/run_certificate.sh \
     circuits/proofbase/sha256_compress_96324_proofbase.xag circuits/sha256_compress_96324.xag "$OUT/twin_cut_96324" "$JOBS"

# 4. compositional proof: proof base == SHA-256 reference.  Quick mode proves all 108 schedule and round lemmas
#    locally and takes the two fused-round tail lemmas from their recorded manifests (the composition check
#    rebuilds those cones and checks their hashes against the manifests); full mode proves all 110.
export VERIFY_TAIL=1 VERIFY_SPLIT=1
W="$OUT/compositional"; mkdir -p "$W"
ONLY=$(python3 -c "print(','.join(['sched:%d'%t for t in range(16,62)]+['round:%d'%t for t in range(62)]))")
comp_quick() {
  cp circuits/proofbase/sha256_compress_96273_proofbase.xag "$W/pb96273.xag"
  rm -rf "$W/rest96273" "$W/merged96273"
  (cd compositional_proof && python3 lemmas.py "$W/pb96273.xag" "$W/rest96273" --jobs "$JOBS" --only "$ONLY" | tail -1 \
   && python3 merge_runs.py "$W/merged96273" "$W/rest96273" ../manifests/tail_lemmas_96273_proofbase \
   && python3 check_compose.py "$W/pb96273.xag" "$W/merged96273/lemmas.json" "$W/merged96273")
}
step compositional_96273_108of110 '"ok": true' comp_quick
# 4b. negative control for the compositional prover: one AND input of the round-10 cone complemented;
#     the round-10 lemma must be reported not_equivalent
comp_neg() {
  python3 tools/mutate_and.py circuits/proofbase/sha256_compress_96273_proofbase.xag "$W/mut96273.xag" 10 || return 1
  rm -rf "$W/mut"
  (cd compositional_proof && python3 lemmas.py "$W/mut96273.xag" "$W/mut" --jobs 1 --only round:10) \
    | tee "$W/mut.log"
  grep -q "^round:10 not_equivalent" "$W/mut.log" && echo "COMPOSITIONAL NEGATIVE CONTROL OK"
}
step compositional_negative_control "COMPOSITIONAL NEGATIVE CONTROL OK" comp_neg
unset VERIFY_TAIL VERIFY_SPLIT

# 5. the reference side of every SAT lemma (329 lemmas: both proof bases and the 100,658 circuit) is node for node
#    the fragment proved equal to the Cryptol specification in fips_spec_proof/ (no SAT, no SAW needed)
step ref_link_329 "ALL REF LINKS OK: 329 lemmas" python3 fips_spec_proof/check_ref_link.py all

# 6. negative controls quoted in the paper: 40 mutated twin-cut lemmas refuted; 21 mutants rejected by both GF(2) checkers
step negative_controls "NEGATIVE CONTROLS OK" sh tools/negative_controls.sh

# 7. variable-CV SHA-256: regenerate the 104,023-gate generator output, and check the shipped 103,669-gate circuit
#    (the same circuit after the joint sigma0/sigma1 linear-layer pass) against it with both GF(2) checkers
varcv() {
  V="$OUT/varcv"; mkdir -p "$V"
  (cd generator && python3 gen_sha2.py "$V/sha256_var_104023.xag" 32 0 > "$V/gen.log" 2>&1 && python3 count.py "$V/sha256_var_104023.xag") || return 1
  # each check's own exit status counts (no pipe hides it); all three must succeed
  python3 generator/simsha2.py "$V/sha256_var_104023.xag" 32 0 10000 > "$V/sim.log" 2>&1; rs=$?
  tail -1 "$V/sim.log"
  python3 rewrite_checks/lincheck.py "$V/sha256_var_104023.xag" circuits/ports/sha256_var_compress_103669.xag > "$V/lincheck.log" 2>&1; rl=$?
  echo "lincheck: $(tail -1 "$V/lincheck.log")"
  (cd rewrite_checks/independent && python3 check2.py "$V/sha256_var_104023.xag" ../../circuits/ports/sha256_var_compress_103669.xag) > "$V/check2.log" 2>&1; rc2=$?
  echo "check2: $(tail -1 "$V/check2.log")"
  if [ "$rs" -ne 0 ]; then echo "simulation FAILED (exit $rs)"; return 1; fi
  if [ "$rl" -ne 0 ] || [ "$rc2" -ne 0 ]; then echo "GF(2) check FAILED (exit $rl / $rc2)"; return 1; fi
  case "$(tail -1 "$V/lincheck.log")" in PASS*) ;; *) return 1;; esac
  case "$(tail -1 "$V/check2.log")" in PASS*) ;; *) return 1;; esac
  echo "VARCV OK"
}
step variable_cv "VARCV OK" varcv

# 8. comparators that can be regenerated offline
step textbook_baseline '"total": 111903' python3 comparators/textbook.py "$OUT/textbook.json"
emp() { python3 comparators/emp_sha256.py fix "$OUT/emp_fix.xag" && python3 comparators/emp_sha256.py var "$OUT/emp_var.xag"; }
step emp_transliteration "fix {'AND': 21863, 'XOR': 89859, 'total': 111722}" emp
if [ -n "${PUBLIC_RAW:-}" ]; then step public_baselines "PUBLIC OK" python3 comparators/score_public.py "$PUBLIC_RAW"; fi

if [ "$MODE" = full ]; then
  G="$OUT/rebuild"; mkdir -p "$G" "$G/cache"
  gen_all() {
    J=$(python3 -c "print(','.join(map(str,range(1,62))))")
    cd generator
    ( GEN_FUSEF=1 GEN_JS=joint_sigma_slp_105.txt MADD2_REG=0 GEN_CACHE="$G/cache" GEN_LAMBDA=8 \
        python3 generate_sha256.py "$G/a_B.xag" $J > "$G/a.log" 2>&1 ) &
    ( GEN_FUSEF=1 GEN_JS=joint_sigma_slp_105.txt MADD2_REG=0 GEN_CACHE="$G/cache" GEN_LAMBDA=8 GEN_MOVE0=1 GEN_J12=1,2 \
        python3 generate_sha256.py "$G/b_B.xag" $J > "$G/b.log" 2>&1 ) &
    ( JOINT_LAMBDA=8 python3 generate_sha256_joint.py "$G/c_raw.xag" final+wfuse+sig105 > "$G/c.log" 2>&1 ) &
    wait
    python3 peep.py "$G/a_B.xag" "$G/a_Bp.xag" && python3 peep.py "$G/b_B.xag" "$G/b_Bp.xag" && python3 peep.py "$G/c_raw.xag" "$G/c.xag"
    cd ..
    cmp "$G/a_Bp.xag" circuits/proofbase/sha256_compress_96273_proofbase.xag && \
    cmp "$G/b_Bp.xag" circuits/proofbase/sha256_compress_96324_proofbase.xag && \
    cmp "$G/c.xag" circuits/sha256_compress_100658.xag && echo "GENERATORS BYTE-IDENTICAL"
  }
  step rebuild_generators "GENERATORS BYTE-IDENTICAL" gen_all
  walks() {
    python3 rewrite_checks/rewrite_check.py "$G/a_Bp.xag" - "$G/a_Bpr.xag" && \
    sh residual_walk/walk.sh "$G/a_Bpr.xag" "$G/a_Bp.xag" "$G/a_pass1.xag" && \
    sh residual_walk/walk.sh "$G/a_pass1.xag" "$G/a_Bp.xag" "$G/a_pass2.xag" 100000000 8 && \
    python3 rewrite_checks/rewrite_check.py "$G/b_Bp.xag" - "$G/b_Bpr.xag" && \
    sh residual_walk/walk.sh "$G/b_Bpr.xag" "$G/b_Bp.xag" "$G/b_pass1.xag" && \
    cmp "$G/a_pass2.xag" circuits/sha256_compress_96273.xag && cmp "$G/b_pass1.xag" circuits/sha256_compress_96324.xag \
    && echo "WALKS BYTE-IDENTICAL"
  }
  step rebuild_walks "WALKS BYTE-IDENTICAL" walks
  comp_full() {  # $1 = circuit file, $2 = tag, $3 = split (1/0)
    cp "$1" "$W/$2.xag"; rm -rf "$W/full_$2"
    (cd compositional_proof && VERIFY_TAIL=1 VERIFY_SPLIT=$3 python3 lemmas.py "$W/$2.xag" "$W/full_$2" --jobs "$JOBS" --timeout 10800 | tail -1 \
     && VERIFY_TAIL=1 VERIFY_SPLIT=$3 python3 check_compose.py "$W/$2.xag" "$W/full_$2/lemmas.json" "$W/full_$2")
  }
  step compositional_96273_full '"ok": true' comp_full circuits/proofbase/sha256_compress_96273_proofbase.xag f96273 1
  step compositional_96324_full '"ok": true' comp_full circuits/proofbase/sha256_compress_96324_proofbase.xag f96324 1
  step compositional_100658_full '"ok": true' comp_full circuits/sha256_compress_100658.xag f100658 0
  ports() {
    cd generator
    python3 gen_sha1.py "$G/sha1.xag" pair+fuse && cmp "$G/sha1.xag" ../circuits/ports/sha1_compress_49851.xag && \
    python3 simsha1.py "$G/sha1.xag" 20000 && \
    python3 gen_sha2.py "$G/sha512_var.xag" 64 0 cshift+neg && python3 simsha2.py "$G/sha512_var.xag" 64 0 5000 && \
    python3 gen_sha2.py "$G/sha512_fix.xag" 64 1 cshift+neg && python3 simsha2.py "$G/sha512_fix.xag" 64 1 5000 && \
    python3 count.py "$G/sha512_var.xag" && python3 count.py "$G/sha512_fix.xag" && echo "PORTS OK"
  }
  step ports_sha1_sha512 "PORTS OK" ports
fi
echo; echo "summary ($MODE mode):"; cat "$SUMMARY"
grep -q FAIL "$SUMMARY" && exit 1 || exit 0
