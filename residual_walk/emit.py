"""emit.py UNIVERSE.xag STATE OUT.xag : rebuild an XAG from a steiner state (v c d per alive XOR).
Every XOR is checked with exact affine forms over the universe's AND basis (vec[v] == vec[c]^vec[d]
up to the constant, which becomes a complemented edge); AND gates are copied.  Raises on any mismatch
or cycle.  The result is then checked independently by lincheck.py."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lin import load, forms, X
sys.setrecursionlimit(1000000)
x = load(sys.argv[1]); vec, _ = forms(x); base = x.n_inputs + 1
defn = {}
for l in open(sys.argv[2]):
    v, c, d = map(int, l.split()); defn[v] = (c, d)
y = X.XAG(x.n_inputs)
lit = {0: 0}
for i in range(1, base): lit[i] = 2 * i
def build(root):
    st = [(root, 0)]; onpath = set()
    while st:
        v, ph = st.pop()
        if v in lit: continue
        op, a, b = x.nodes[v - base]
        ops = defn.get(v, (a >> 1, b >> 1)) if op == "AND" else defn[v]
        if ph == 0:
            assert v not in onpath, "cycle"
            onpath.add(v); st.append((v, 1))
            for u in ops:
                if u not in lit: st.append((u, 0))
        else:
            onpath.discard(v)
            if op == "AND":
                F = vec[a >> 1] ^ (a & 1); G = vec[b >> 1] ^ (b & 1)
                c, d = ops
                def coords(P):
                    for a1, b1 in ((1, 0), (0, 1), (1, 1)):
                        r = P ^ (F if a1 else 0) ^ (G if b1 else 0)
                        if r in (0, 1): return a1, b1, r
                    raise AssertionError(f"AND {v}: operand outside plane")
                a1, b1, c1 = coords(vec[c]); a2, b2, c2 = coords(vec[d])
                assert (a1 & b2) ^ (a2 & b1) == 1, f"AND {v}: dependent operands"
                done = False
                for al in (0, 1):
                    for be in (0, 1):
                        k1, k2 = c1 ^ al, c2 ^ be
                        if ((a1 & a2) ^ (a1 & k2) ^ (a2 & k1)) == 0 and ((b1 & b2) ^ (b1 & k2) ^ (b2 & k1)) == 0 and not done:
                            lit[v] = y.AND(lit[c] ^ al, lit[d] ^ be) ^ (k1 & k2); done = True
                assert done
            else:
                c, d = ops
                k = vec[v] ^ vec[c] ^ vec[d]
                assert k in (0, 1), f"form mismatch at {v}"
                lit[v] = y.XOR(lit[c], lit[d]) ^ k
for o in x.outputs: build(o >> 1)
y.outputs = [lit[o >> 1] ^ (o & 1) for o in x.outputs]
open(sys.argv[3], "w").write(X.write_xag(y))
print(y.counts())
