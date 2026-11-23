#!/usr/bin/env python3
"""abc_control.py RAWDIR [--abc PATH] [--deepsyn] -- the ABC control comparators of the paper's Table 2.

RAWDIR must hold both raw files (NIST sha256-input512.slp and Campanelli et al. sha256Unpadded.txt); a missing
or altered file, an ABC failure, a failed simulation or an unexpected result gives exit status 1.

Recipe (historical run of 2026-11-18: ABC commit c2c7046d built with the default Makefile and
ABC_USE_NO_READLINE=1 in an Ubuntu 24.04 container on x86-64 (AMD EPYC 9B14, Google Cloud c3d), one run per
pipeline, single-threaded):
  1. the public file (NIST sha256-input512.slp or Campanelli et al. sha256Unpadded.txt, hash-checked against
     ../RAW_SHA256SUMS) is put in SPEC order (interface found by simulation) and canonicalised with
     verifier/cost.py: NIST 113,386 gates, Campanelli et al. 113,797;
  2. it is written as BLIF (xagio.xag_to_blif) and run through ABC:
       map only     read_blif IN; strash; read_genlib xag.genlib; map -a; write_blif OUT
       compress2rs  read_blif IN; strash; <compress2rs script, once>; read_genlib xag.genlib; map -a; write_blif OUT
       deepsyn      read_blif IN; strash; &get -n; &deepsyn -T 600 -S 1; &put; read_genlib xag.genlib; map -a; ...
     (xag.genlib: AND2, XOR2, XNOR2 area 1; INV, BUF and constants area 0);
  3. the mapped BLIF is read back (xagio.blif_to_xag), canonicalised and scored with verifier/cost.py, and
     simulated against the FIPS 180-4 reference.
Expected canonical results (AND, XOR, total, canonical sha256) are in EXPECTED.  The Campanelli et al. map-only
row depends on the ABC build: the same commit built on macOS (Apple clang) gives a different, also correct,
mapping (ALSO_OBSERVED), while a Linux arm64 build with gcc 13.3 reproduces the historical result;
either recorded observation is accepted for that row, and the output says which one matched.  Every circuit also matched the
reference on 100,004 simulated inputs in our run (this script simulates 2,000).  These are smoke-tested
comparators, not proved circuits.  Selection rule of the paper: the smallest canonical total among all
measured fixed-IV comparators, including these controls, is the headline comparator (here 109,780).
"""
import argparse
import hashlib
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "verifier"))
sys.path.insert(0, HERE)
import cost  # noqa: E402
import sim  # noqa: E402
from xag import load  # noqa: E402
from xagio import xag_to_blif, blif_to_xag  # noqa: E402

C2RS = ("balance -l; resub -K 6 -l; rewrite -l; resub -K 6 -N 2 -l; refactor -l; resub -K 8 -l; balance -l; "
        "resub -K 8 -N 2 -l; rewrite -l; resub -K 10 -l; rewrite -z -l; resub -K 10 -N 2 -l; balance -l; "
        "resub -K 12 -l; refactor -z -l; resub -K 12 -N 2 -l; rewrite -z -l; balance -l")
PIPES = {"mapA": "", "compress2rs": C2RS + "; ", "deepsyn": "&get -n; &deepsyn -T 600 -S 1; &put; "}
INPUTS = {"nist": "sha256-input512.slp", "campanelli": "sha256Unpadded.txt"}
EXPECTED = {
    ("nist", "mapA"): (22296, 89152, 111448, "17902aee486deaf82f3827ce37e7e4aa45ed19a6c5d830f7f91095d612c2d012"),
    ("nist", "compress2rs"): (62580, 48348, 110928, "380d52dba08153a00e8da7794360a9df838211061c6e461edce3ecc1e67a8665"),
    ("nist", "deepsyn"): (61978, 48336, 110314, "bdba2a97d5547aa1d7324c22d1a717c2a2bb26b47a7affc850401943eaff4f8a"),
    ("campanelli", "mapA"): (22272, 89261, 111533, "1fa2ec35964b19effc1802cfe5638390701ec881ff81eb2e8482e9477e840af7"),
    ("campanelli", "compress2rs"): (61356, 48424, 109780, "b68730ca7fe99898c1a2a05b2bba5628a2324e9809d0dd9764a0a3287f2d9c72"),
    ("campanelli", "deepsyn"): (62086, 48325, 110411, "f957fb56a5334193e958f02a238d20bacf5b507602b0f5f8ff4e1854317e3041"),
}
# Other recorded observations, accepted for the named row only (AND, XOR, total, canonical sha256 or its prefix).
ALSO_OBSERVED = {
    ("campanelli", "mapA"): [("macOS arm64 build (Apple clang 21) of the same ABC commit", 22272, 89224, 111496,
                            "8fd4a25dac7d37858ed50b77d1d8c99c502d8e3656ba998f273f783bd74869a5")],
}


def spec_ordered(path):
    x = load(path).to_xag()
    hits = sim.discover(x, functions=["primary"])
    if not hits:
        raise SystemExit("%s: interface not found" % path)
    h = hits[0]
    return cost.canonicalize(sim.normalize(x, h["in_map"], h["out_map"], 512))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rawdir")
    ap.add_argument("--abc", default="abc")
    ap.add_argument("--deepsyn", action="store_true", help="also run the 10-minute deepsyn pipelines")
    a = ap.parse_args()
    sums = {}
    for ln in open(os.path.join(HERE, "..", "RAW_SHA256SUMS")):
        h, p = ln.split()
        sums.setdefault(os.path.basename(p), set()).add(h)
    tmp = tempfile.mkdtemp()
    genlib = os.path.join(HERE, "xag.genlib")
    ok = True
    for tag, fn in INPUTS.items():
        p = os.path.join(a.rawdir, fn)
        if not os.path.isfile(p):   # the default replay needs both raw input files
            print("%-10s %s missing: FAILED" % (tag, fn)); ok = False; continue
        if hashlib.sha256(open(p, "rb").read()).hexdigest() not in sums.get(fn, ()):
            print("%-10s %s does not match RAW_SHA256SUMS" % (tag, fn)); ok = False; continue
        y = spec_ordered(p)
        r0 = cost.report_xag(y)
        print("%-10s input canonical AND %d XOR %d total %d" % (tag, r0["AND"], r0["XOR"], r0["total"]))
        blif = os.path.join(tmp, tag + ".blif")
        open(blif, "w").write(xag_to_blif(y))
        for pipe, script in PIPES.items():
            if pipe == "deepsyn" and not a.deepsyn:
                continue
            out = os.path.join(tmp, "%s_%s.blif" % (tag, pipe))
            cmd = "read_blif %s; strash; %sread_genlib %s; map -a; write_blif %s" % (blif, script, genlib, out)
            rc = subprocess.run([a.abc, "-c", cmd], capture_output=True, text=True).returncode
            if rc != 0 or not os.path.exists(out):
                print("    %-12s ABC failed (exit %d)" % (pipe, rc)); ok = False; continue
            z = cost.canonicalize(blif_to_xag(open(out).read()))
            r = cost.report_xag(z)
            v = sim.verify(z, "primary", list(range(512)), [(i, 0) for i in range(256)], n=2000)
            got = (r["AND"], r["XOR"], r["total"], cost.xag_hash(z))
            note = None
            if got == EXPECTED[(tag, pipe)]:
                note = "OK" + (" (matches the historical Linux x86-64 run)" if (tag, pipe) in ALSO_OBSERVED else "")
            for what, *e in ALSO_OBSERVED.get((tag, pipe), []):
                if got[:3] == tuple(e[:3]) and got[3].startswith(e[3]):
                    note = "OK (matches the %s; platform-dependent row)" % what
            good = note is not None and v["status"] == "pass"
            ok &= good
            print("    %-12s AND %d XOR %d total %d sha256 %s simulation %s  %s"
                  % (pipe, got[0], got[1], got[2], got[3][:16], v["status"],
                     note if good else "DIFFERENT from the recorded runs %s" % (EXPECTED[(tag, pipe)][:3],)))
    print("ABC CONTROL %s" % ("OK" if ok else "DIFFERENT OR FAILED"))
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
