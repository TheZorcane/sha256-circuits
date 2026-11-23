"""joint-sum planner generator: constant-redistribution generator gen4 structure (chpair, final, wfuse) with the three round sums T1, a', e' replaced
by ONE joint two-output column compression (jmadd.build2), which contains column-compression generator's "inline", constant-redistribution generator's constant
shift (solved exactly inside the DP) and the classical shared-T1 round as special cases.
Word-level cuts W_t, a_t, e_t are unchanged (T1 no longer exists as a word in general).

usage: [JOINT_LAMBDA=x] [JOINT_MODE=0|1] generate_sha256_joint.py OUT.xag [variant]     variant tokens: final, wfuse, neg, sig105
"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen3
from gen3 import XAG, FB, g2, sref, cost, chpair, write_xag
import jmadd

M32 = (1 << 32) - 1
# linear-layer optimiser's joint sigma0/sigma1 SLP (results/sigma_ctx_105.txt, copied unmodified): lines "mask p q", mask = p ^ q
SLP = [tuple(map(int, l.split())) for l in open(os.path.join(HERE, "joint_sigma_slp_105.txt"))]


def build(variant="final+wfuse"):
    x = XAG(512); b = FB(x)
    PC1, PC2 = {}, {}
    S1 = lambda rows, K=0: jmadd.build1(x, rows, K, cache=PC1)
    W = g2.inputs_words(x, 16)
    Wrows = {}
    sig = {}
    def sigmas(t):
        """(sigma0(W_t), sigma1(W_t)); with 'sig105' the 35 words that need both use linear-layer optimiser's joint 105-XOR SLP."""
        if t not in sig:
            if "sig105" in variant and 14 <= t <= 48:
                node = {1 << i: W[t][i] for i in range(32)}
                for g, p, q in SLP: node[g] = b.xor(node[p], node[q])
                tap = lambda k, rots, sh: sum(1 << ((k + r) % 32) for r in rots) | ((1 << (k + sh)) if k + sh < 32 else 0)
                sig[t] = ([node[tap(k, (7, 18), 3)] for k in range(32)], [node[tap(k, (17, 19), 10)] for k in range(32)])
            else:
                sig[t] = (b.s0(W[t]) if t <= 48 else None, b.s1(W[t]) if t >= 14 else None)
        return sig[t]
    for t in range(16, 64):
        Wrows[t] = [sigmas(t-2)[1], W[t-7], sigmas(t-15)[0], W[t-16]]
        W.append(None if ("wfuse" in variant and t >= 62) else S1(Wrows[t]))
    st = [g2.const_word(v) for v in sref.IV]
    log = []
    for t in range(64):
        a, bb, c, d, e, f, gg, h = st
        wr = [W[t]] if W[t] is not None else Wrows[t]
        t1rows, t1pairs = [b.S1(e)] + wr, [(h, chpair(x, e, f, gg, h))]
        last = (t == 63 and "final" in variant)
        ka = sref.IV[0] if last else 0
        ke = sref.IV[4] if last else 0
        s0, mj = b.S0(a), b.maj(a, bb, c)
        info = {}
        n0 = len(x.nodes); a0 = sum(1 for nd in x.nodes if nd[0] == "AND")
        neg = False
        if "neg" in variant:      # negate the whole joint heap: complement every row and both outputs
            nS, nA, nE = len(t1rows) + 2 * len(t1pairs), 2, 1
            nrows = [[l ^ 1 for l in r] for r in t1rows]
            npairs = [([u ^ 1 for u in U], Wd) for (U, Wd) in t1pairs]
            ns0, nmj, nd = [l ^ 1 for l in s0], [l ^ 1 for l in mj], [l ^ 1 for l in d]
            kS, kA, kE = (nS - 1 - sref.K[t]) & M32, (nA - ka) & M32, (nE - ke) & M32
            c_pos = jmadd.estimate2(t1rows, t1pairs, [s0, mj], [d], sref.K[t], ka, ke, cache=PC2)
            c_neg = jmadd.estimate2(nrows, npairs, [ns0, nmj], [nd], kS, kA, kE, cache=PC2)
            neg = c_neg < c_pos
        if neg:
            an, en = jmadd.build2(x, nrows, npairs, [ns0, nmj], [nd], kS, kA, kE, cache=PC2, info=info)
            an, en = [l ^ 1 for l in an], [l ^ 1 for l in en]
        else:
            an, en = jmadd.build2(x, t1rows, t1pairs, [s0, mj], [d], sref.K[t], ka, ke, cache=PC2, info=info)
        a1 = sum(1 for nd_ in x.nodes if nd_[0] == "AND")
        log.append((t, info["cost"], hex(info["cS"]), "neg" if neg else "pos", "AND", a1 - a0, "gates", len(x.nodes) - n0,
                    "partial_cols", [i for i, r in enumerate(info["rem"]) if r is not None and r > 1]))
        st = [an, a, bb, c, en, e, f, gg]
    if "final" in variant:
        out = [st[0]] + [S1([st[k]], sref.IV[k]) for k in (1, 2, 3)] + [st[4]] + \
              [S1([st[k]], sref.IV[k]) for k in (5, 6, 7)]
    else:
        out = [S1([st[k]], sref.IV[k]) for k in range(8)]
    x.outputs = [w[31 - j] for w in out for j in range(32)]
    return x, log


if __name__ == "__main__":
    variant = sys.argv[2] if len(sys.argv) > 2 else "final+wfuse"
    x, log = build(variant)
    for l in log: print(*l, file=sys.stderr)
    y = cost.canonicalize(x)
    r = cost.report_xag(y)
    r["lambda"] = jmadd.LAM; r["joint"] = jmadd.JOINT
    print(json.dumps(r))
    if len(sys.argv) > 1: open(sys.argv[1], "w").write(write_xag(y))
