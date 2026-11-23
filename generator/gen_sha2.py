"""port generator generic SHA-2 compression generator (SHA-256 / SHA-512, fixed IV or variable chaining value).
Adaptation of constant-redistribution generator gen4.py (column-compression generator gen3 + constant redistribution) to word size WS and variable CV.
usage: gen_sha2.py OUT.xag WS FIXED(0|1) [variant]     variant tokens: cshift neg final wfuse
Inputs (SPEC order): [CV words (variable only)] + 16 block words, each word MSB first.
Outputs: 8 words H_out, each MSB first."""
import sys, os, json, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen3
from gen3 import XAG, X, A, cost, PC, write_xag
import madd
from sha2ref import Sha2


def build(ws, fixed, variant="cshift+neg+final+wfuse", iters=int(os.environ.get("GEN_ITERS", "60"))):
    R = Sha2(ws); M = R.M; NR = R.R
    rng = random.Random(20261118 + ws)
    nin = 16 * ws + (0 if fixed else 8 * ws)
    x = XAG(nin)
    base = 0 if fixed else 8 * ws

    def words(first, n):
        return [[x.input(first + ws * k + ws - 1 - j) for j in range(ws)] for k in range(n)]
    CV = None if fixed else words(0, 8)
    W = words(base, 16)

    def S(rows, K=0, pairs=()): return madd.build_sum(x, rows, K & M, n=ws, plan_cache=PC, pairs=pairs)
    def E(rows, K=0, pairs=()): return madd.estimate(rows, K & M, n=ws, pairs=pairs, plan_cache=PC)
    def rotr(w, n): return [w[(i + n) % ws] for i in range(ws)]
    def shr(w, n): return [w[i + n] if i + n < ws else 0 for i in range(ws)]
    def xorw(*ws_):
        r = ws_[0]
        for w in ws_[1:]: r = [X(x, p, q) for p, q in zip(r, w)]
        return r
    def sig(w, t): return xorw(rotr(w, t[0]), rotr(w, t[1]), rotr(w, t[2]))
    def ssg(w, t): return xorw(rotr(w, t[0]), rotr(w, t[1]), shr(w, t[2]))
    def maj(a, b, c): return [X(x, b[i], A(x, X(x, a[i], b[i]), X(x, b[i], c[i]))) for i in range(ws)]

    def negate(rows, K, pairs):
        n = len(rows) + 2 * len(pairs)
        return [[l ^ 1 for l in r] for r in rows], (n - 1 - K) & M, [([u ^ 1 for u in U], Wd) for (U, Wd) in pairs]

    def En(rows, K, pairs=(), neg=False):
        if neg: rows, K, pairs = negate(rows, K, pairs)
        return E(rows, K, pairs)

    def Sn(rows, K, pairs=(), neg=False):
        if neg:
            r2, K2, p2 = negate(rows, K, pairs)
            return [l ^ 1 for l in S(r2, K2, p2)]
        return S(rows, K, pairs)

    def choose(t1rows, t1pairs, Kt, s0, mj, d, use_neg, ka=0, ke=0):
        V = [2] * ws
        NEG = list(itertools.product((False, True), repeat=3)) if use_neg else [(False,) * 3]
        def rc(C, negs):
            return (En(t1rows, Kt - C, t1pairs, negs[0]) + En([V, s0, mj], C + ka, (), negs[2]) +
                    En([d, V], C + ke, (), negs[1]))
        def best_neg(C): return min((rc(C, n), n) for n in NEG)
        bc = 0; bv, bn = best_neg(0); base = bv
        for it in range(iters):
            C = rng.getrandbits(ws) if it % 2 else bc ^ (1 << rng.randrange(ws))
            v, n = best_neg(C)
            if v < bv: bv, bc, bn = v, C, n
        imp = True
        while imp:
            imp = False
            for i in range(ws):
                for C in (bc ^ (1 << i), (bc + (1 << i)) & M, (bc - (1 << i)) & M):
                    v, n = best_neg(C)
                    if v < bv: bv, bc, bn, imp = v, C, n, True
        return base, bv, bc, bn

    def chpair(x_, e, f, g, h):
        return [X(x, X(x, g[i], h[i]), A(x, e[i], X(x, f[i], g[i]))) for i in range(ws)]

    log = []
    Wrows = {}
    for t in range(16, NR):
        Wrows[t] = [ssg(W[t-2], R.ss1), W[t-7], ssg(W[t-15], R.ss0), W[t-16]]
        W.append(None if ("wfuse" in variant and t >= NR - 2) else S(Wrows[t]))
    if fixed: st = [[(v >> i) & 1 for i in range(ws)] for v in R.IV]
    else: st = list(CV)
    for t in range(NR):
        a, bb, c, d, e, f, gg, h = st
        wr = [W[t]] if W[t] is not None else Wrows[t]
        t1rows, t1pairs = [sig(e, R.bs1)] + wr, [(h, chpair(x, e, f, gg, h))]
        last = (t == NR - 1 and "final" in variant)
        ka = R.IV[0] if (last and fixed) else 0
        ke = R.IV[4] if (last and fixed) else 0
        xa = [CV[0]] if (last and not fixed) else []
        xe = [CV[4]] if (last and not fixed) else []
        Kt = R.K[t]
        s0, mj = sig(a, R.bs0), maj(a, bb, c)
        V = [2] * ws
        sep = E(t1rows, Kt, t1pairs) + E([V, s0, mj] + xa, ka) + E([d, V] + xe, ke)
        inl = E(t1rows + [s0, mj] + xa, Kt + ka, t1pairs) + E(t1rows + [d] + xe, Kt + ke, t1pairs)
        if inl < sep:
            an = S(t1rows + [s0, mj] + xa, Kt + ka, t1pairs)
            en = S(t1rows + [d] + xe, Kt + ke, t1pairs)
            log.append((t, "inline", sep, inl))
        elif last and not fixed:
            t1 = S(t1rows, Kt, t1pairs)
            an = S([t1, s0, mj] + xa); en = S([d, t1] + xe)
            log.append((t, "last-sep", sep, inl))
        else:
            C, negs = 0, (False,) * 3
            if "cshift" in variant or "neg" in variant:
                base_, bv, bc, bn = choose(t1rows, t1pairs, Kt, s0, mj, d, "neg" in variant, ka, ke)
                if "cshift" not in variant: bc = 0
                C, negs = bc, bn
                log.append((t, base_, bv, hex(bc), bn))
            t1 = Sn(t1rows, Kt - C, t1pairs, negs[0])
            an = Sn([t1, s0, mj], C + ka, (), negs[2])
            en = Sn([d, t1], C + ke, (), negs[1])
        st = [an, a, bb, c, en, e, f, gg]
    out = []
    for k in range(8):
        if "final" in variant and k in (0, 4):
            out.append(st[k]); continue
        out.append(S([st[k]], R.IV[k]) if fixed else S([st[k], CV[k]]))
    x.outputs = [w[ws - 1 - j] for w in out for j in range(ws)]
    return x, log


if __name__ == "__main__":
    ws = int(sys.argv[2]); fixed = int(sys.argv[3])
    variant = sys.argv[4] if len(sys.argv) > 4 else "cshift+neg+final+wfuse"
    x, log = build(ws, fixed, variant)
    for l in log[-3:]: print(*l, file=sys.stderr)
    y = cost.canonicalize(x)
    print(json.dumps(cost.report_xag(y)))
    open(sys.argv[1], "w").write(write_xag(y))
