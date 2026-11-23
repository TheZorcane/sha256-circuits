"""Simulate a SHA-2 compression .xag (SPEC order) against sha2ref.
usage: simsha2.py F.xag WS FIXED(0|1) [n] [seed]
Exit status 0 only if every vector matches."""
import sys, os, json, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
import sim
from xag import read_xag
from sha2ref import Sha2
x = read_xag(open(sys.argv[1]).read())
ws = int(sys.argv[2]); fixed = int(sys.argv[3])
n = int(sys.argv[4]) if len(sys.argv) > 4 else 100000
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 12345
R = Sha2(ws); M = R.M
r = random.Random(seed)
edge = [0, M, 1, 1 << (ws - 1), 0x5555555555555555 & M, 0xaaaaaaaaaaaaaaaa & M, M >> 1, (1 << (ws - 1)) + 1]
def rw(): return r.choice(edge) if r.random() < .15 else r.getrandbits(ws)
nin = (16 if fixed else 24)
vs = [[rw() for _ in range(nin)] for _ in range(n)]
vs[0] = [0] * nin; vs[1] = [M] * nin
# known answer: padded "abc" block (with IV as CV for variable)
blk = R.pad(b"abc")[0]
vs[2] = blk if fixed else list(R.IV) + blk
t = time.time()
width = ws * nin
lns_ = sim.vectors_to_lanes([R.words_to_int(v) for v in vs], width)
out = sim.lanes_to_vectors(sim.simulate(x, lns_, n), n)
bad = 0
for o, v in zip(out, vs):
    H, B = (list(R.IV), v) if fixed else (v[:8], v[8:])
    if o != R.words_to_int(R.compress(H, B)): bad += 1
print(json.dumps({"file": sys.argv[1], "ws": ws, "fixed_iv": bool(fixed), "vectors": n, "seed": seed,
                  "mismatches": bad, "status": "pass" if bad == 0 else "FAIL", "seconds": round(time.time() - t, 1)}))
sys.exit(0 if bad == 0 else 1)   # nonzero process status on any mismatch
