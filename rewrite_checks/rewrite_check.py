"""rewrite_check.py B.xag CAND.xag|- [OUT.xag] : third, purely structural link between the proof variant B and the candidate.

Rebuild B gate by gate, applying only the XOR identity   p ^ (q ^ p) = q   (with edge complements) whenever
an XOR gate has an operand that is itself an XOR gate containing the other operand; everything else is copied.
Then canonicalize with verifier/cost.py and compare with CAND.  Soundness: every rewritten node computes the
same function as the node of B it replaces (induction over the file order), so the rebuilt circuit is
equivalent to B; if its canonical text equals CAND's, CAND is that circuit.
Independent of lincheck (no linear algebra, no AND matching)."""
import sys, os, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
from xag import XAG, read_xag, write_xag
import cost

b = read_xag(open(sys.argv[1]).read())
y = XAG(b.n_inputs)
base = b.n_inputs + 1
lit = list(range(0, 2 * base, 2))        # lit[v] = literal in y of node v of b (inputs map to themselves)
kind = {}                                # y-node -> (a, b) literals if it is an XOR gate of y
nrew = 0
for j, (op, a, c) in enumerate(b.nodes):
    la = lit[a >> 1] ^ (a & 1); lc = lit[c >> 1] ^ (c & 1)
    if op == "AND":
        r = y.AND(la, lc)
    else:
        r = None
        for p, q in ((la, lc), (lc, la)):
            d = kind.get(q >> 1)
            if d is None: continue
            u, w = d
            if u >> 1 == p >> 1: r = w ^ (u & 1) ^ (p & 1) ^ (q & 1); break
            if w >> 1 == p >> 1: r = u ^ (w & 1) ^ (p & 1) ^ (q & 1); break
        if r is None:
            r = y.XOR(la, lc)
            if r >> 1 >= base and (r >> 1) not in kind: kind[r >> 1] = (la, lc) if not (r & 1) else (la ^ 1, lc)
        else:
            nrew += 1
    lit.append(r)
y.outputs = [lit[o >> 1] ^ (o & 1) for o in b.outputs]
yc = cost.canonicalize(y)
t = write_xag(yc)
cand = open(sys.argv[2]).read() if sys.argv[2] != "-" else ""
h1 = hashlib.sha256(t.encode()).hexdigest(); h2 = hashlib.sha256(cand.encode()).hexdigest()
r = cost.report_xag(yc)
print("rewrites applied:", nrew, "| rebuilt:", r["AND"], r["XOR"], r["total"], h1[:16], "| candidate:", h2[:16])
print("PASS: rebuilt circuit is byte-identical to the candidate" if t == cand else "DIFFERENT (no claim)")
if len(sys.argv) > 3: open(sys.argv[3], "w").write(t)
