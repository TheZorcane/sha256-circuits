"""Simulate a SHA-1 compression .xag (512 -> 160, SPEC order) against sha1ref; optional interface discovery of the
old-Bristol sha-1.txt.  usage: simsha1.py F.xag [n] [seed]   |   simsha1.py --bristol F.txt"""
import sys, os, json, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
import sim, xag
import sha1ref as R


def vecs(n, seed):
    r = random.Random(seed)
    edge = [0, R.M, 1, 1 << 31, 0x55555555, 0xaaaaaaaa, 0x7fffffff, 0x80000001]
    vs = [[r.choice(edge) if r.random() < .15 else r.getrandbits(32) for _ in range(16)] for _ in range(n)]
    vs[0] = [0] * 16; vs[1] = [R.M] * 16; vs[2] = R.pad(b"abc")[0]
    return vs


if sys.argv[1] == "--bristol":
    rc = xag.load(sys.argv[2]); x = rc.to_xag()
    N = 200; vs = vecs(N, 3)
    refo = [R.words_to_bits(R.compress(list(R.IV), v)) for v in vs]
    def rw(b): return [b[32 * (i // 32) + 31 - i % 32] for i in range(len(b))]
    found = []
    for irev in ("none", "rev_all", "rev_word"):
        cols = []
        for v in vs:
            b = R.words_to_bits(v)
            cols.append(b if irev == "none" else b[::-1] if irev == "rev_all" else rw(b))
        o = sim.simulate(x, [sum(cols[k][c] << k for k in range(N)) for c in range(512)], N)
        for orev in ("none", "rev_all", "rev_word"):
            ok = True
            for k in range(N):
                ob = [(o[c] >> k) & 1 for c in range(160)]
                ob = ob if orev == "none" else ob[::-1] if orev == "rev_all" else rw(ob)
                if ob != refo[k]: ok = False; break
            if ok: found.append((irev, orev))
    print("bristol sha-1 interface (input,output) =", found)
    sys.exit(0)

x = xag.read_xag(open(sys.argv[1]).read())
n = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 12345
vs = vecs(n, seed)
t = time.time()
lns_ = sim.vectors_to_lanes([R.words_to_int(v) for v in vs], 512)
out = sim.lanes_to_vectors(sim.simulate(x, lns_, n), n)
bad = sum(1 for o, v in zip(out, vs) if o != R.words_to_int(R.compress(list(R.IV), v)))
print(json.dumps({"file": sys.argv[1], "vectors": n, "seed": seed, "mismatches": bad,
                  "status": "pass" if bad == 0 else "FAIL", "seconds": round(time.time() - t, 1)}))
sys.exit(0 if bad == 0 else 1)   # nonzero process status on any mismatch
