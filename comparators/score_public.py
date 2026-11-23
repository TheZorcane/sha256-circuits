#!/usr/bin/env python3
"""score_public.py RAWDIR -- score the public SHA-256 circuits under the paper's cost model.

RAWDIR holds the files downloaded by fetch_public_baselines.sh.  Each file's sha256 must match RAW_SHA256SUMS
(the exact files scored for the paper).  For each file we print the file's own gate tally, the direct NOT-free
AND+XOR count, and the canonical count (verifier/cost.py canonicalize: NOT free, XNOR = XOR, constants folded,
structural hashing, dead gates removed).  Expected canonical totals: NIST fixed-IV 113,386 (22,296 AND),
Goldfeder unpadded 113,797 (22,272 AND), NIST variable-CV 116,239, Bristol Fashion sha256 132,319.
For the two variable-CV files it also derives a fixed-IV circuit (interface discovered by simulation, chaining
value tied to the standard IV, canonicalised, simulated on 2,000 vectors): expected NIST 113,386 (22,296 AND, the
same netlist as the fixed-IV file) and Bristol Fashion 129,422 (22,356 AND)."""
import os, sys, json, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "verifier"))
import cost
import sim
import sha256_ref as ref
from xag import load

DERIVED = {"sha256-iv256+input512.slp": (22296, 91090, 113386), "sha256.txt": (22356, 107066, 129422)}


def derived_fixed_iv(rc):
    """Tie the chaining-value inputs of a variable-CV circuit to the standard IV (derived, not published)."""
    x = rc.to_xag()
    hits = sim.discover(x, functions=["secondary"])
    if len(hits) != 1:
        return None, "no unique interface"
    h = hits[0]
    y = cost.canonicalize(sim.normalize(x, h["in_map"], h["out_map"], sim.FUNCTIONS["secondary"][0]))
    in_map = [("const", b) for b in ref.words_to_bits(list(ref.IV))] + list(range(512))
    z = cost.canonicalize(sim.normalize(y, in_map, [(i, 0) for i in range(256)], 512))
    v = sim.verify(z, "primary", list(range(512)), [(i, 0) for i in range(256)], n=2000)
    return z, v["status"]


FILES = {
    "sha256-input512.slp": ("NIST SLP, fixed IV (512 -> 256)", (22296, 91090, 113386), "slp"),
    "sha256Unpadded.txt": ("Goldfeder / CCS 2017, fixed IV", (22272, 91525, 113797), "old"),
    "sha256-iv256+input512.slp": ("NIST SLP, variable CV (768 -> 256)", (22573, 93666, 116239), "slp"),
    "sha256.txt": ("Bristol Fashion sha256, variable CV", (22573, 109746, 132319), "fashion"),
}


def main(rawdir):
    sums = {}
    for ln in open(os.path.join(HERE, "RAW_SHA256SUMS")):
        h, p = ln.split()
        sums.setdefault(os.path.basename(p), set()).add(h)
    ok = True; seen = 0
    for fn, (desc, exp, fmt) in FILES.items():
        p = os.path.join(rawdir, fn)
        if not os.path.exists(p):
            print("%-30s missing (skipped)" % fn); continue
        seen += 1
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        if h not in sums.get(fn, ()):
            print("%-30s sha256 %s does NOT match RAW_SHA256SUMS" % (fn, h)); ok = False; continue
        rep = cost.report(load(p, fmt))
        c = rep["canonical_xag"]; d = rep["direct_xag"]
        good = (c["AND"], c["XOR"], c["total"]) == exp
        ok &= good
        print("%-30s %s\n    file tally %s\n    direct NOT-free AND %d XOR %d total %d | canonical AND %d XOR %d total %d  %s"
              % (fn, desc, json.dumps(rep["raw_counts_by_type"]), d["AND"], d["XOR"], d["total"], c["AND"], c["XOR"], c["total"],
                 "OK" if good else "MISMATCH (expected %s)" % (exp,)))
        if fn in DERIVED:
            z, st = derived_fixed_iv(load(p, fmt))
            if z is None:
                print("    derived fixed IV: %s" % st); ok = False; continue
            cz = cost.report_xag(z)
            good = (cz["AND"], cz["XOR"], cz["total"]) == DERIVED[fn] and st == "pass"
            ok &= good
            print("    derived fixed IV (CV tied to the IV, canonical): AND %d XOR %d total %d, simulation %s  %s"
                  % (cz["AND"], cz["XOR"], cz["total"], st, "OK" if good else "MISMATCH (expected %s)" % (DERIVED[fn],)))
    print("PUBLIC %s (%d files scored)" % ("OK" if ok and seen else "FAIL" if not ok else "NONE", seen))
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1]) else 1)
