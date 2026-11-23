"""twin-cut certificate lemma library, version 2: twin cut points.

Cut points are *twin pairs* (u, t): u a vnew node, t a vref node with u == t ^ pol as functions of the 512 primary
inputs.  Pair search uses exact GF(2) forms over (inputs, vref AND outputs) -- search only; every pair is re-proved
by an LRAT-certified SAT lemma:

  lemma(u, t, pol):  free variables = vref node ids (the leaf set);
     A = gate of u over its operands, where every strictly-lower vnew node s that has a twin T[s] is replaced
         by free var T[s] (^ polarity), other vnew nodes expanded gate by gate down to such nodes;
     B = gate of t over its operands, every vref node in TW (the set of twin targets) replaced by its free var,
         others (incl. all vref ANDs not in TW) expanded gate by gate.
     claim: A ^ pol == B  for all values of the free variables.
Induction in vnew topological order: A-leaf nodes were proved earlier, so evaluating at true values gives
value_new(u) == value_ref(t) ^ pol.
Var ids are XAG var ids (0 const, 1..n inputs, n+1+k node k).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
from xag import XAG, AND, XOR, read_xag


def load(path):
    return read_xag(open(path).read())


def all_forms(X, n, repl=None, R=None, fr=None):
    """affine forms (python ints; bit0 const, bit 1+i input i, bit 1+n+j  vref AND j) per var id.
    repl given -> X is vnew, ANDs replaced by plane relation to vref (R,fr)."""
    rb = n + 1
    f = [0] + [1 << (1 + i) for i in range(n)]
    lf = lambda l: f[l >> 1] ^ (l & 1)
    j = 0
    jof = {}
    for idx, (op, a, b) in enumerate(X.nodes):
        if op == XOR:
            f.append(lf(a) ^ lf(b))
        elif repl is None:
            f.append(1 << (1 + n + j)); jof[rb + idx] = j; j += 1
        else:
            jv, cf, cg, c = repl[rb + idx]
            F = fr[R.nodes[jv - rb][1] >> 1] ^ (R.nodes[jv - rb][1] & 1)
            G = fr[R.nodes[jv - rb][2] >> 1] ^ (R.nodes[jv - rb][2] & 1)
            f.append((1 << (1 + n + RJ[jv])) ^ (F if cf else 0) ^ (G if cg else 0) ^ c)
    return f, jof


RJ = {}


def match(R, N):
    """plane matching, same algebra as lincheck: returns M: N AND var -> (R AND var, cf, cg, c)"""
    n = R.n_inputs
    rb, nb = n + 1, N.n_inputs + 1
    fr, jof = all_forms(R, n)
    RJ.update(jof)
    key = {}
    for idx, (op, a, b) in enumerate(R.nodes):
        if op == AND:
            f = fr[a >> 1] ^ (a & 1); g = fr[b >> 1] ^ (b & 1)
            f0, g0 = f & ~1, g & ~1
            if f0 and g0 and f0 != g0:
                key.setdefault(tuple(sorted((f0, g0, f0 ^ g0))), (rb + idx, f, g))
    form = [0] + [1 << (1 + i) for i in range(n)]
    lf = lambda l: form[l >> 1] ^ (l & 1)
    M = {}
    for idx, (op, a, b) in enumerate(N.nodes):
        if op == XOR:
            form.append(lf(a) ^ lf(b)); continue
        p, q = lf(a), lf(b)
        p0, q0 = p & ~1, q & ~1
        hit = key.get(tuple(sorted((p0, q0, p0 ^ q0)))) if p0 and q0 and p0 != q0 else None
        if hit is None:
            raise SystemExit(f"no plane for N node {idx}")
        jv, f, g = hit

        def co(P):
            for a1, b1 in ((1, 0), (0, 1), (1, 1)):
                r = P ^ (f if a1 else 0) ^ (g if b1 else 0)
                if r in (0, 1):
                    return a1, b1, r
            raise SystemExit("coords")
        a1, b1, k1 = co(p); a2, b2, k2 = co(q)
        assert (a1 & b2) ^ (a2 & b1) == 1
        cf = (a1 & a2) ^ (a1 & k2) ^ (a2 & k1)
        cg = (b1 & b2) ^ (b1 & k2) ^ (b2 & k1)
        M[nb + idx] = (jv, cf, cg, k1 & k2)
        form.append((1 << (1 + n + RJ[jv])) ^ (f if cf else 0) ^ (g if cg else 0) ^ (k1 & k2))
    assert [fr[o >> 1] ^ (o & 1) for o in R.outputs] == [lf(o) for o in N.outputs], "algebraic output mismatch"
    return M, fr, form


def twins(R, N, fr, fn):
    """T: vnew var (input/xor/and) -> (vref var t, pol).  Outputs list: (u, t_out, pol)."""
    first = {}
    for v, f in enumerate(fr):
        if v == 0:
            continue
        first.setdefault(f & ~1, v)
    T = {}
    nb = N.n_inputs + 1
    for v in range(1, len(fn)):
        t = first.get(fn[v] & ~1)
        if t is not None:
            T[v] = (t, (fn[v] ^ fr[t]) & 1)
    outs = []
    for o, (lo, ro) in enumerate(zip(N.outputs, R.outputs)):
        u, t = lo >> 1, ro >> 1
        if u == 0 or t == 0:
            outs.append((o, u, t, None)); continue
        assert (fn[u] & ~1) == (fr[t] & ~1), f"output {o} form differs"
        outs.append((o, u, t, (fn[u] ^ fr[t]) & 1))
    return T, outs




# ---------------------------------------------------------------------------------------------------------
# lemma construction with depth parameter d.
#   A side: vnew node u expanded (its gate), every strictly lower vnew node s with a twin T[s]=(w,pol) is a leaf that
#           stands for vref node w ^ pol; w is replaced by its vref expansion (below).
#   vref expansion set X(d): all vref non-input nodes at distance < d from the roots {t} U {w : A-leaf twins};
#           nodes of X are built as vref gates in both circuits (where reachable), all other referenced vref nodes
#           are the shared free variables.
# ---------------------------------------------------------------------------------------------------------
def item_struct(R, N, T, u, t, d):
    rb, nb = R.n_inputs + 1, N.n_inputs + 1
    an = set(); wa = set()
    stack = [(u, True)]
    while stack:
        v, top = stack.pop()
        if v == 0:
            continue
        if not top and v in T:
            wa.add(T[v][0]); continue
        if v in an:
            continue
        assert v >= nb, "input as lemma top"
        an.add(v)
        _, a, b = N.nodes[v - nb]
        stack.append((a >> 1, False)); stack.append((b >> 1, False))
    roots = set(wa); roots.add(t)
    dist = {}
    front_ = [r for r in roots if r >= rb]
    for r in front_:
        dist[r] = 0
    level = 0
    X = set()
    while front_ and level < d:
        nxt = []
        for v in front_:
            X.add(v)
            _, a, b = R.nodes[v - rb]
            for l in (a, b):
                w = l >> 1
                if w >= rb and w not in dist and w not in X:
                    dist[w] = level + 1; nxt.append(w)
        front_ = nxt; level += 1
    return sorted(an), wa, X


def build(R, N, T, items):
    """items: list of (u, t, pol, d).  Per item private gate copies.  Returns (A, B, free)."""
    rb, nb = R.n_inputs + 1, N.n_inputs + 1
    st = [item_struct(R, N, T, u, t, d) for u, t, pol, d in items]
    # free = vref vars (inputs or non-X nodes) referenced
    free = set()
    for (u, t, pol, d), (an, wa, X) in zip(items, st):
        refs = set(wa); refs.add(t)
        for v in sorted(X):
            _, a, b = R.nodes[v - rb]
            refs.add(a >> 1); refs.add(b >> 1)
        free |= {v for v in refs if v != 0 and (v < rb or v not in X)}
    free = sorted(free)
    fidx = {v: i for i, v in enumerate(free)}
    A = XAG(len(free)); B = XAG(len(free))
    oA = []; oB = []

    def expand(Xc, X, roots):
        """build X-gates of vref reachable from roots in circuit Xc; returns map var->lit for built + free"""
        need = set(); stack = [r for r in roots if r in X]
        while stack:
            v = stack.pop()
            if v in need:
                continue
            need.add(v)
            _, a, b = R.nodes[v - rb]
            for l in (a, b):
                if (l >> 1) in X:
                    stack.append(l >> 1)
        m = {0: 0}
        for v in sorted(need):
            op, a, b = R.nodes[v - rb]
            la = m[a >> 1] ^ (a & 1) if (a >> 1) in m else 2 * (fidx[a >> 1] + 1) ^ (a & 1)
            lb = m[b >> 1] ^ (b & 1) if (b >> 1) in m else 2 * (fidx[b >> 1] + 1) ^ (b & 1)
            m[v] = A_or_B_gate(Xc, op, la, lb)
        return m

    def A_or_B_gate(Xc, op, la, lb):
        return Xc.AND(la, lb) if op == AND else Xc.XOR(la, lb)

    for (u, t, pol, d), (an, wa, X) in zip(items, st):
        xa = expand(A, X, wa)
        xb = expand(B, X, [t])
        vl = lambda m, w: m[w] if w in m else 2 * (fidx[w] + 1)
        am = {0: 0}
        for v in an:
            op, a, b = N.nodes[v - nb]

            def lit(l):
                w = l >> 1
                if w in am:
                    return am[w] ^ (l & 1)
                tt, pl = T[w]
                return vl(xa, tt) ^ pl ^ (l & 1)
            am[v] = A_or_B_gate(A, op, lit(a), lit(b))
        oA.append(am[u] ^ pol)
        oB.append(vl(xb, t))
    A.outputs = oA; B.outputs = oB
    A.validate(); B.validate()
    return A, B, free


def sim_equal(A, B, rounds=4, seed=1):
    """fast random simulation of A vs B on shared free vars (256 vectors per round); True if all outputs equal."""
    import random
    rnd = random.Random(seed)
    W = 256; mask = (1 << W) - 1
    for _ in range(rounds):
        vals = [0] + [rnd.getrandbits(W) for _ in range(A.n_inputs)]
        res = []
        for X in (A, B):
            v = list(vals)
            for op, a, b in X.nodes:
                x = v[a >> 1] ^ (mask if a & 1 else 0)
                y = v[b >> 1] ^ (mask if b & 1 else 0)
                v.append(x & y if op == AND else x ^ y)
            res.append([v[o >> 1] ^ (mask if o & 1 else 0) for o in X.outputs])
        if res[0] != res[1]:
            return [i for i in range(len(res[0])) if res[0][i] != res[1][i]]
    return []
