"""augment.py BASE.xag QUADS.txt OUT.xag : append, as dead candidate nodes, the pair sums of every
non-trivial weight-4 dependency (Q lines) and of P lines (explicit pairs)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lin import load, forms, X
x = load(sys.argv[1]); vec, _ = forms(x)
have = {v & ~1 for v in vec}
new = {}
for l in open(sys.argv[2]):
    t = l.split()
    if not t: continue
    if t[0] == 'Q':
        a, b, c, d = map(int, t[1:5]); prs = [(a, b), (a, c), (a, d)]
    elif t[0] == 'P':
        prs = [(int(t[1]), int(t[2]))]
    else: continue
    for p, q in prs:
        f = (vec[p] ^ vec[q]) & ~1
        if f and f not in have and f not in new: new[f] = (p, q)
for f, (p, q) in new.items(): x.XOR(2 * p, 2 * q)
if len(sys.argv) > 4:
    # AND planes: AND(f,g) may take any two of {f, g, f^g} as operands (polarities fix the function).
    # Append the missing thirds as dead candidates and write "P N f g h" lines.
    vec, _ = forms(x); base = x.n_inputs + 1
    idx = {}
    for v in range(1, len(vec)): idx.setdefault(vec[v] & ~1, v)
    n0 = len(x.nodes); out = []
    for j in range(n0):
        op, a, b = x.nodes[j]
        if op != "AND": continue
        f, g = a >> 1, b >> 1
        hf = (vec[f] ^ vec[g]) & ~1
        if hf == 0 or f == 0 or g == 0: continue
        if hf not in idx:
            idx[hf] = (x.XOR(2 * f, 2 * g)) >> 1
        out.append(f"P {base + j} {f} {g} {idx[hf]}\n")
    open(sys.argv[4], "w").writelines(out)
    print("planes", len(out), "thirds added", len(x.nodes) - n0)
open(sys.argv[3], "w").write(X.write_xag(x))
print("added", len(new), "nodes; total", len(x.nodes))
