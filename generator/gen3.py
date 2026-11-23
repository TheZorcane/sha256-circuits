"""column-compression generator (fixed IV): base generator structure (Maj XOR sharing, Ch = g^(e&(f^g)))
with every multi-operand sum built by madd.build_sum (MDFA column compression,
DP-planned).  usage: gen3.py OUT.xag"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
from xag import XAG, write_xag
import cost
import sha256_ref as sref
import gen as g2
import madd

PC = {}


def S(x, rows, K=0, pairs=()):
    return madd.build_sum(x, rows, K, plan_cache=PC, pairs=pairs)


def X(x, a, b):
    if a in (0, 1): return b ^ a
    if b in (0, 1): return a ^ b
    if a >> 1 == b >> 1: return (a ^ b) & 1
    return x.XOR(a, b)


def A(x, a, b):
    if a in (0, 1): return b if a else 0
    if b in (0, 1): return a if b else 0
    if a == b: return a
    if a == b ^ 1: return 0
    return x.AND(a, b)


class FB(g2.B):
    """base generator's word builder with constant folding (so madd sees constant bits as constants)."""
    def xor(s, a, b): return X(s.x, a, b)
    def and_(s, a, b): return A(s.x, a, b)


def chpair(x, e, f, g, h):
    """pair word (h, h^Ch): h^Ch = (g^h) ^ (e & (f^g)); g^h is round t-1's f^g (strash-shared)."""
    return [X(x, X(x, g[i], h[i]), A(x, e[i], X(x, f[i], g[i]))) for i in range(32)]


def E(rows, K=0, pairs=()):
    return madd.estimate(rows, K, pairs=pairs, plan_cache=PC)


def build(variant="v1"):
    x = XAG(512); b = FB(x)
    W = g2.inputs_words(x, 16)
    Wrows = {}
    for t in range(16, 64):
        rows = [b.s1(W[t-2]), W[t-7], b.s0(W[t-15]), W[t-16]]
        Wrows[t] = rows
        W.append(S(x, rows) if not ("wfuse" in variant and t >= 62) else None)
    st = [g2.const_word(v) for v in sref.IV]
    for t in range(64):
        a, bb, c, d, e, f, gg, h = st
        wr = Wrows[t] if W[t] is None else [W[t]]
        if "chpair" in variant:
            t1rows, t1pairs = [b.S1(e)] + wr, [(h, chpair(x, e, f, gg, h))]
        else:
            t1rows, t1pairs = [h, b.S1(e), b.ch(e, f, gg)] + wr, []
        Kt = sref.K[t]
        s0, mj = b.S0(a), b.maj(a, bb, c)
        last = (t == 63 and "final" in variant)
        Ka = sref.IV[0] if last else 0
        Ke = sref.IV[4] if last else 0
        inline = False
        if "inline" in variant:
            sep = E(t1rows, Kt, t1pairs) + E([[2] * 32, s0, mj], Ka) + E([d, [2] * 32], Ke)
            inl = E(t1rows + [s0, mj], Kt + Ka, t1pairs) + E(t1rows + [d], Kt + Ke, t1pairs)
            inline = inl < sep
            if t < 6 or inline: print("round", t, sep, inl, file=sys.stderr)
        if inline:
            an = S(x, t1rows + [s0, mj], Kt + Ka, t1pairs)
            en = S(x, t1rows + [d], Kt + Ke, t1pairs)
        else:
            t1 = S(x, t1rows, Kt, t1pairs)
            an = S(x, [t1, s0, mj], Ka)
            en = S(x, [d, t1], Ke)
        st = [an, a, bb, c, en, e, f, gg]
    if "final" in variant:
        out = [st[0]] + [S(x, [st[k]], sref.IV[k]) for k in (1, 2, 3)] + [st[4]] + \
              [S(x, [st[k]], sref.IV[k]) for k in (5, 6, 7)]
    else:
        out = [S(x, [st[k]], sref.IV[k]) for k in range(8)]
    x.outputs = [w[31 - j] for w in out for j in range(32)]
    return x


if __name__ == "__main__":
    variant = sys.argv[2] if len(sys.argv) > 2 else "v1"
    x = build(variant)
    y = cost.canonicalize(x)
    r = cost.report_xag(y)
    print(json.dumps(r))
    if len(sys.argv) > 1:
        open(sys.argv[1], "w").write(write_xag(y))
