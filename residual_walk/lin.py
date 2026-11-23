"""linear-layer optimiser: linear-layer extraction.  Every node of an XAG is an affine form over the basis
{1, primary inputs, AND outputs}.  vec[var] is a python int: bit 0 = constant, bit 1+i = PI i,
bit 1+n_in+j = j-th AND gate (in file order)."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
import xag as X


def load(path):
    return X.read_xag(open(path).read())


def forms(x):
    n = x.n_inputs
    base = n + 1
    vec = [0] * (base + len(x.nodes))
    for i in range(n):
        vec[i + 1] = 1 << (1 + i)
    andidx = {}
    k = 0
    for j, (op, a, b) in enumerate(x.nodes):
        v = base + j
        if op == "AND":
            andidx[v] = k
            vec[v] = 1 << (1 + n + k)
            k += 1
        else:
            vec[v] = vec[a >> 1] ^ vec[b >> 1] ^ ((a ^ b) & 1)
    return vec, andidx


if __name__ == "__main__":
    import collections
    x = load(sys.argv[1])
    vec, andidx = forms(x)
    n = x.n_inputs; base = n + 1
    print(x.counts())
    tg = {}
    for j, (op, a, b) in enumerate(x.nodes):
        if op == "AND":
            for l in (a, b):
                tg.setdefault(vec[l >> 1] & ~1, []).append(j)
    for o in x.outputs:
        tg.setdefault(vec[o >> 1] & ~1, []).append(-1)
    triv = sum(1 for f in tg if bin(f).count("1") <= 1)
    print("distinct target forms", len(tg), "trivial (basis)", triv)
    xorforms = {}
    for j, (op, a, b) in enumerate(x.nodes):
        if op == "XOR":
            xorforms.setdefault(vec[base + j] & ~1, []).append(j)
    print("xor nodes", sum(len(v) for v in xorforms.values()), "distinct forms", len(xorforms))
    print("xor nodes that are targets", sum(1 for f in xorforms if f in tg))
    ws = sorted(bin(f).count("1") for f in tg)
    print("target weight quantiles", [ws[int(q * (len(ws) - 1))] for q in (0, .1, .25, .5, .75, .9, 1)])
    fo = collections.Counter()
    for j, (op, a, b) in enumerate(x.nodes):
        fo[a >> 1] += 1; fo[b >> 1] += 1
    for o in x.outputs: fo[o >> 1] += 1
    inter = [j for j, (op, a, b) in enumerate(x.nodes) if op == "XOR" and (vec[base + j] & ~1) not in tg]
    print("intermediate xor nodes", len(inter), "fanout hist",
          sorted(collections.Counter(fo[base + j] for j in inter).items()))
    tx = [j for j, (op, a, b) in enumerate(x.nodes) if op == "XOR" and (vec[base + j] & ~1) in tg]
    print("target xor nodes fanout hist", sorted(collections.Counter(fo[base + j] for j in tx).items()))
