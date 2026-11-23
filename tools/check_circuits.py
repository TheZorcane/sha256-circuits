#!/usr/bin/env python3
"""check_circuits.py -- recount, simulate and export the circuits of the paper.

  check_circuits.py count      canonical AND/XOR counts, depth and canonical sha256 of every .xag in circuits/;
                               asserts the numbers quoted in the paper
  check_circuits.py sim [N]    random-vector simulation (default N = 100,000 random + 4 known vectors) of the SHA-256
                               circuits (fixed IV, and the variable-CV circuit in ports/) against the independent FIPS 180-4 reference in verifier/sha256_ref.py
                               (a smoke test, not a proof)
  check_circuits.py bristol    write Bristol Fashion files (verifier/xag.py write_bristol_fashion) into
                               circuits/bristol/ and re-read each one with the Bristol parser: the re-read
                               circuit must canonicalize to the same counts and pass simulation
  check_circuits.py bristol-check   only the re-read + simulation of the shipped Bristol files
Run with PYTHONPATH=../verifier (reproduce.sh does this)."""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "verifier"))
import cost, sim
from xag import read_xag, write_bristol_fashion, load
sys.path.insert(0, os.path.join(HERE, "..", "comparators"))
from rescore import not_tally

CI = os.path.join(HERE, "..", "circuits")
EXPECT = {  # canonical AND, XOR, total
    "sha256_compress_96273.xag": (21783, 74490, 96273),
    "sha256_compress_96324.xag": (21746, 74578, 96324),
    "sha256_compress_100658.xag": (21745, 78913, 100658),
    "proofbase/sha256_compress_96273_proofbase.xag": None,
    "proofbase/sha256_compress_96324_proofbase.xag": None,
    "proofbase/walk_pass1.xag": (21783, 74491, 96274),
    "ports/sha1_compress_49851.xag": (10552, 39299, 49851),
    "ports/sha256_var_compress_103669.xag": (22293, 81376, 103669),
}
SHA256 = [k for k in EXPECT if not k.startswith("ports/")]
NOT_TALLY = {  # NOT-counting tally comparable to NIST's (XOR may be emitted as XNOR): total, NOT
    "sha256_compress_96273.xag": (99035, 2762),
    "sha256_compress_96324.xag": (99054, 2730),
}


def count():
    ok = True
    for name, exp in EXPECT.items():
        x = read_xag(open(os.path.join(CI, name)).read())
        c = cost.canonicalize(x)
        r = cost.report_xag(c)
        same = cost.report_xag(x)["sha256"] == r["sha256"]
        good = exp is None or (r["AND"], r["XOR"], r["total"]) == exp
        ok &= good
        print("%-48s AND %6d XOR %6d total %6d depth %5d AND-depth %4d canonical=%s %s"
              % (name, r["AND"], r["XOR"], r["total"], r["depth"], r["and_depth"], "yes" if same else "no",
                 "OK" if good else "MISMATCH"))
    for name, exp in NOT_TALLY.items():
        t = not_tally(read_xag(open(os.path.join(CI, name)).read()))
        good = (t["tally_total"], t["NOT"]) == exp
        ok &= good
        print("%-48s NOT-counting tally %d (%d NOT) %s" % (name, t["tally_total"], t["NOT"], "OK" if good else "MISMATCH"))
    print("COUNTS %s" % ("OK" if ok else "FAIL"))
    return ok


def simulate(n):
    ok = True
    for name in SHA256:
        x = read_xag(open(os.path.join(CI, name)).read())
        v = sim.verify(x, "primary", list(range(x.n_inputs)), [(i, 0) for i in range(256)], n=n, seed=20261119)
        print("%-48s %s vectors=%d mismatches=%d" % (name, v["status"], v["vectors_tested"], v["mismatches"]))
        ok &= v["mismatches"] == 0
    name = "ports/sha256_var_compress_103669.xag"     # variable chaining value: 768 inputs (CV words, then block)
    x = read_xag(open(os.path.join(CI, name)).read())
    v = sim.verify(x, "secondary", list(range(x.n_inputs)), [(i, 0) for i in range(256)], n=n, seed=20261119)
    print("%-48s %s vectors=%d mismatches=%d" % (name, v["status"], v["vectors_tested"], v["mismatches"]))
    ok &= v["mismatches"] == 0
    print("SIM %s" % ("OK" if ok else "FAIL"))
    return ok


def bristol(write):
    ok = True
    BD = os.path.join(CI, "bristol")
    os.makedirs(BD, exist_ok=True)
    for name in ("sha256_compress_96273.xag", "sha256_compress_96324.xag", "sha256_compress_100658.xag"):
        x = read_xag(open(os.path.join(CI, name)).read())
        bp = os.path.join(BD, name.replace(".xag", ".bristol"))
        if write:
            # one input group of 512 bits (block), one output group of 256 bits; SPEC bit order
            open(bp, "w").write(write_bristol_fashion(x, [512], [256]))
        y = cost.canonicalize(load(bp, "fashion").to_xag()) if True else None
        r0, r1 = x.counts(), y.counts()
        v = sim.verify(y, "primary", list(range(y.n_inputs)), [(i, 0) for i in range(256)], n=2000, seed=7)
        good = (r0["AND"], r0["XOR"]) == (r1["AND"], r1["XOR"]) and v["mismatches"] == 0
        ok &= good
        print("%-48s re-read AND %d XOR %d, sim %s mismatches %d %s" % (os.path.basename(bp), r1["AND"], r1["XOR"], v["status"], v["mismatches"], "OK" if good else "FAIL"))
    print("BRISTOL %s" % ("OK" if ok else "FAIL"))
    return ok


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "count":
        r = count()
    elif m == "sim":
        r = simulate(int(sys.argv[2]) if len(sys.argv) > 2 else 100000)
    elif m == "bristol":
        r = bristol(True)
    elif m == "bristol-check":
        r = bristol(False)
    else:
        raise SystemExit(__doc__)
    sys.exit(0 if r else 1)
