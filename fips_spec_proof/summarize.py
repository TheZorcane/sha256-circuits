"""Check and summarize the logs of fips_spec_proof/run.sh.  usage: summarize.py LOGDIR

Exit status 0 only if every required result is present and correct:
  env.txt            cryptol-specs commit is exactly SPEC_COMMIT
  gen.log            120 / 109 / 110 fragments (plain / tail / split), inline identity, validated mutant
  frags_{S,I}_{plain,tail,split}.log   (all six)  SAW exit 0, exactly N distinct "PROVED F<style>_<var>_*"
                     theorems (N = 120 / 109 / 110), the "ALL N FRAGMENT THEOREMS" line, and the same
                     fragment names in styles S and I
  compose_{plain,tail,split}.log       (all three)  SAW exit 0, "sanity OK", PROVED W_I (composition),
                     PROVED B_S and B_S_0 .. B_S_447 (449 bridge theorems), PROVED B_S_ALL_LENGTHS:
                     exactly 451 distinct PROVED lines
  negative.log       SAW exit 0 and "NEGATIVE CONTROL OK" (the mutant was rejected with a counterexample)
  ref_link.log       329 per-lemma OK lines, the three per-circuit "REF LINK OK" lines, and
                     "ALL REF LINKS OK: 329 lemmas" as the last line
A missing log, a missing or duplicated theorem, a wrong count, a nonzero SAW exit, a timeout or a failure
line in any required log makes the exit status 1.  Other *.log files in LOGDIR are listed and ignored.
"""
import json
import os
import re
import sys

SPEC_COMMIT = "8b74f5fc2b77fc7bbb587ba2217990b1963e9d6c"
N_FRAG = {"plain": 120, "tail": 109, "split": 110}
REF_LINES = ["REF LINK OK: 110 of 110 lemmas (sha256_compress_96273_proofbase.xag, variant split)",
             "REF LINK OK: 110 of 110 lemmas (sha256_compress_96324_proofbase.xag, variant split)",
             "REF LINK OK: 109 of 109 lemmas (sha256_compress_100658.xag, variant tail)"]
BAD = re.compile(r"TIMEOUT|MEMORY LIMIT|FAILED|unsolved subgoal|^Invalid|^prove: |^Stack trace|error",
                 re.M | re.I)

failures = []


def fail(msg):
    failures.append(msg)
    print("    FAIL: " + msg)


def read(L, name):
    p = os.path.join(L, name)
    if not os.path.isfile(p):
        fail(f"{name} missing")
        return None
    return open(p, errors="replace").read()


def saw_log(L, name, allow_bad=False):
    """common checks for a run_saw.py log; returns (text, list of PROVED names) or (None, [])"""
    txt = read(L, name)
    if txt is None:
        return None, []
    m = re.search(r"^# run_saw: exit=(\d+) wall=([0-9.]+)s", txt, re.M)
    proved = re.findall(r"^PROVED (\S+)\s*$", txt, re.M)
    times = [float(t) for t in re.findall(r"^Time: ([0-9.]+)s", txt, re.M)]
    print(f"{name:20s} exit={m.group(1) if m else None} wall={m.group(2) if m else None}s "
          f"proved={len(proved)} solver_time_sum={sum(times):.1f}s")
    if not m:
        fail(f"{name}: no run_saw exit line (run interrupted?)")
    elif m.group(1) != "0":
        fail(f"{name}: SAW exit status {m.group(1)}")
    if not allow_bad:
        bad = [ln for ln in txt.splitlines() if BAD.search(ln)]
        if bad:
            fail(f"{name}: failure line: {bad[0].strip()[:120]}")
    if len(set(proved)) != len(proved):
        fail(f"{name}: duplicated PROVED lines")
    return txt, proved


def main(L):
    # environment: exact spec commit
    env = read(L, "env.txt")
    if env is not None:
        m = re.search(r"^cryptol-specs commit: (\S+)", env, re.M)
        got = m.group(1) if m else None
        print(f"env.txt              cryptol-specs commit {got}")
        if got != SPEC_COMMIT:
            fail(f"env.txt: cryptol-specs commit {got}, expected {SPEC_COMMIT}")
    # generator
    g = read(L, "gen.log")
    if g is not None:
        try:
            gi = json.loads(g)
            ok = ({k: gi.get(k) for k in N_FRAG} == N_FRAG
                  and gi["full_reference"]["identical_to_inlined_fragments_plain"] is True
                  and gi["mutant_round_10"]["aig_validated_vs_sim"] is True)
        except (ValueError, KeyError, TypeError):
            ok = False
        print("gen.log              " + ("fragment counts 120/109/110, inline identity, mutant validated"
                                         if ok else "BAD"))
        if not ok:
            fail("gen.log: fragment counts, inline identity or mutant validation not as expected")
    # six fragment logs
    names = {}
    for var, n in N_FRAG.items():
        for st in ("S", "I"):
            name = f"frags_{st}_{var}.log"
            txt, proved = saw_log(L, name)
            if txt is None:
                continue
            pre = f"F{st}_{var}_"
            own = {p[len(pre):] for p in proved if p.startswith(pre)}
            if len(proved) != n or len(own) != n:
                fail(f"{name}: {len(own)} distinct {pre}* theorems proved, expected {n}")
            if f"ALL {n} FRAGMENT THEOREMS F{st} PROVED FOR VARIANT {var}" not in txt:
                fail(f"{name}: final 'ALL {n} FRAGMENT THEOREMS' line missing")
            names[(st, var)] = own
        if ("S", var) in names and ("I", var) in names and names[("S", var)] != names[("I", var)]:
            fail(f"frags_*_{var}.log: styles S and I prove different fragment sets")
    # three composition logs
    want = {"W_I", "B_S", "B_S_ALL_LENGTHS"} | {f"B_S_{k}" for k in range(448)}
    for var in N_FRAG:
        name = f"compose_{var}.log"
        txt, proved = saw_log(L, name)
        if txt is None:
            continue
        if f"== variant {var}" not in txt:
            fail(f"{name}: not a log of variant {var}")
        if not re.search(r"^sanity OK\s*$", txt, re.M):
            fail(f"{name}: 'sanity OK' missing")
        missing = want - set(proved)
        if missing or len(proved) != len(want):
            fail(f"{name}: {len(set(proved) & want)} of {len(want)} required theorems proved"
                 + (f" (missing e.g. {sorted(missing)[0]})" if missing else f" ({len(proved)} PROVED lines)"))
    # negative control (its log legitimately contains the word 'counterexample')
    txt, _ = saw_log(L, "negative.log", allow_bad=True)
    if txt is not None:
        if "NEGATIVE CONTROL FAILED" in txt or not re.search(r"^NEGATIVE CONTROL OK", txt, re.M):
            fail("negative.log: negative control did not reject the mutant")
        else:
            print("    NEGATIVE CONTROL OK")
    # reference link
    r = read(L, "ref_link.log")
    if r is not None:
        lines = [ln.rstrip() for ln in r.splitlines() if ln.strip()]
        per = [ln for ln in lines if re.match(r"^\S+\s+(plain|tail|split)\s+nodes", ln)]
        okper = [ln for ln in per if ln.endswith("  OK")]
        print(f"ref_link.log         {len(okper)} of {len(per)} per-lemma links OK; last line: "
              f"{lines[-1] if lines else '(empty)'}")
        if len(per) != 329 or len(okper) != 329:
            fail(f"ref_link.log: {len(okper)} OK per-lemma lines of {len(per)}, expected 329 of 329")
        for want_ln in REF_LINES:
            if want_ln not in lines:
                fail(f"ref_link.log: missing '{want_ln}'")
        if not lines or lines[-1] != "ALL REF LINKS OK: 329 lemmas":
            fail("ref_link.log: last line is not 'ALL REF LINKS OK: 329 lemmas'")
    required = {"env.txt", "gen.log", "negative.log", "ref_link.log"} | \
               {f"frags_{s}_{v}.log" for s in "SI" for v in N_FRAG} | {f"compose_{v}.log" for v in N_FRAG}
    extra = sorted(f for f in os.listdir(L) if f.endswith(".log") and f not in required) if os.path.isdir(L) else []
    if extra:
        print("not required (ignored): " + ", ".join(extra))
    if failures:
        print(f"SOME REQUIRED CHECK FAILED ({len(failures)} failures)")
        return 1
    print("ALL REQUIRED CHECKS PASSED: 678 fragment theorems, 3 composition and 3 x 449 bridge theorems, "
          "negative control, 329 reference links, cryptol-specs " + SPEC_COMMIT[:12])
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    sys.exit(main(sys.argv[1]))
