"""Simulate a PRIMARY .xag (SPEC order) against verifier/sha256_ref.
usage: simcheck.py F.xag [n] [seed]"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
import sim
from xag import read_xag
x = read_xag(open(sys.argv[1]).read())
n = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 12345
in_map = list(range(512))
out_map = [(i, 0) for i in range(256)]
v = sim.verify(x, "primary", in_map, out_map, n=n, seed=seed)
print(json.dumps(v))
sys.exit(0 if v["mismatches"] == 0 and v["status"] == "pass" else 1)
