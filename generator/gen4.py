"""constant-redistribution generator: column-compression generator gen3 (chpair+inline) + per-round constant redistribution and sum negation.
  T1' = T1 + K_t - C ;  e' = d + T1' + C ;  a' = T1' + S0 + Maj + C        (C chosen per round)
  any sum may be negated: complement all rows/pairs and the output, constant c -> n-1-c.
Word-level cuts W_t, a_t, e_t are unchanged (only T1 changes value).
usage: gen4.py OUT.xag [variant]   variant tokens: cshift, neg, final, wfuse"""
import sys, os, json, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen3
from gen3 import XAG, FB, g2, sref, cost, chpair, madd, PC, write_xag

M32 = (1 << 32) - 1


def nrows(rows, pairs): return len(rows) + 2 * len(pairs)


def negate(rows, K, pairs):
    n = nrows(rows, pairs)
    r2 = [[l ^ 1 for l in r] for r in rows]
    p2 = [([u ^ 1 for u in U], Wd) for (U, Wd) in pairs]
    return r2, (n - 1 - K) & M32, p2


def E(rows, K, pairs=(), neg=False):
    if neg: rows, K, pairs = negate(rows, K, pairs)
    return madd.estimate(rows, K & M32, pairs=pairs, plan_cache=PC)


def S(x, rows, K, pairs=(), neg=False):
    if neg:
        r2, K2, p2 = negate(rows, K, pairs)
        return [l ^ 1 for l in madd.build_sum(x, r2, K2, plan_cache=PC, pairs=p2)]
    return madd.build_sum(x, rows, K & M32, plan_cache=PC, pairs=pairs)


def round_cost(t1rows, t1pairs, Kt, s0, mj, d, C, negs, ka=0, ke=0):
    V = [2] * 32
    return (E(t1rows, Kt - C, t1pairs, negs[0]) + E([V, s0, mj], C + ka, (), negs[2]) +
            E([d, V], C + ke, (), negs[1]))


def choose(t1rows, t1pairs, Kt, s0, mj, d, use_neg, rng, iters=int(os.environ.get("GEN_ITERS", "60")), ka=0, ke=0):
    NEG = list(itertools.product((False, True), repeat=3)) if use_neg else [(False,) * 3]
    def best_neg(C):
        return min((round_cost(t1rows, t1pairs, Kt, s0, mj, d, C, n, ka, ke), n) for n in NEG)
    bc = 0; bv, bn = best_neg(0); base = bv
    for it in range(iters):
        C = rng.getrandbits(32) if it % 2 else bc ^ (1 << rng.randrange(32))
        v, n = best_neg(C)
        if v < bv: bv, bc, bn = v, C, n
    imp = True
    while imp:
        imp = False
        for i in range(32):
            for C in (bc ^ (1 << i), (bc + (1 << i)) & M32, (bc - (1 << i)) & M32):
                v, n = best_neg(C)
                if v < bv: bv, bc, bn, imp = v, C, n, True
    return base, bv, bc, bn


def build(variant="cshift"):
    rng = random.Random(20261118)
    x = XAG(512); b = FB(x)
    W = g2.inputs_words(x, 16)
    Wrows = {}
    for t in range(16, 64):
        Wrows[t] = [b.s1(W[t-2]), W[t-7], b.s0(W[t-15]), W[t-16]]
        W.append(None if ("wfuse" in variant and t >= 62) else gen3.S(x, Wrows[t]))
    st = [g2.const_word(v) for v in sref.IV]
    log = []
    for t in range(64):
        a, bb, c, d, e, f, gg, h = st
        wr = [W[t]] if W[t] is not None else Wrows[t]     # wfuse: W62/W63 rows straight into T1
        t1rows, t1pairs = [b.S1(e)] + wr, [(h, chpair(x, e, f, gg, h))]
        last = (t == 63 and "final" in variant)          # final: fuse IV0/IV4 feed-forward
        ka = sref.IV[0] if last else 0
        ke = sref.IV[4] if last else 0
        Kt = sref.K[t]
        s0, mj = b.S0(a), b.maj(a, bb, c)
        sep = gen3.E(t1rows, Kt, t1pairs) + gen3.E([[2] * 32, s0, mj], ka) + gen3.E([d, [2] * 32], ke)
        inl = gen3.E(t1rows + [s0, mj], Kt + ka, t1pairs) + gen3.E(t1rows + [d], Kt + ke, t1pairs)
        if inl < sep:     # column-compression generator "inline" (fires in round 0)
            an = gen3.S(x, t1rows + [s0, mj], Kt + ka, t1pairs)
            en = gen3.S(x, t1rows + [d], Kt + ke, t1pairs)
            log.append((t, "inline", sep, inl))
        else:
            C, negs = 0, (False,) * 3
            if "cshift" in variant or "neg" in variant:
                base, bv, bc, bn = choose(t1rows, t1pairs, Kt, s0, mj, d, "neg" in variant, rng, ka=ka, ke=ke)
                if "cshift" not in variant: bc = 0
                C, negs = bc, bn
                log.append((t, base, bv, hex(bc), bn))
            t1 = S(x, t1rows, Kt - C, t1pairs, negs[0])
            an = S(x, [t1, s0, mj], C + ka, (), negs[2])
            en = S(x, [d, t1], C + ke, (), negs[1])
        st = [an, a, bb, c, en, e, f, gg]
    if "final" in variant:
        out = [st[0]] + [gen3.S(x, [st[k]], sref.IV[k]) for k in (1, 2, 3)] + [st[4]] + \
              [gen3.S(x, [st[k]], sref.IV[k]) for k in (5, 6, 7)]
    else:
        out = [gen3.S(x, [st[k]], sref.IV[k]) for k in range(8)]
    x.outputs = [w[31 - j] for w in out for j in range(32)]
    return x, log


if __name__ == "__main__":
    variant = sys.argv[2] if len(sys.argv) > 2 else "cshift+neg"
    x, log = build(variant)
    for l in log: print(*l, file=sys.stderr)
    y = cost.canonicalize(x)
    r = cost.report_xag(y)
    print(json.dumps(r))
    if len(sys.argv) > 1: open(sys.argv[1], "w").write(write_xag(y))
