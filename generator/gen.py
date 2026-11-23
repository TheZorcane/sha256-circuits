"""base generator generator: SHA-256 round logic and full PRIMARY circuit as XAG.
Words are lists of 32 literals, LSB first.  Uses verifier/ read-only."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
from xag import XAG
import cost
import sha256_ref as ref

M = 0xFFFFFFFF
F, T = 0, 1


def const_word(v): return [(v >> i) & 1 for i in range(32)]
def rotr(w, n): return [w[(i + n) % 32] for i in range(32)]
def shr(w, n): return [w[i + n] if i + n < 32 else F for i in range(32)]


class B:
    def __init__(s, x): s.x = x
    def xor(s, a, b): return s.x.XOR(a, b)
    def and_(s, a, b): return s.x.AND(a, b)

    def xorw(s, *ws):
        r = ws[0]
        for w in ws[1:]: r = [s.xor(a, b) for a, b in zip(r, w)]
        return r

    def fa(s, a, b, c):
        cs = [v for v in (a, b, c) if v in (0, 1)]
        if cs:
            vs = [v for v in (a, b, c) if v not in (0, 1)]
            k = 0
            for v in cs: k ^= v
            if len(cs) == 1:
                p, q = vs
                sm = s.xor(p, q) ^ k
                cy = s.and_(p ^ 1, q ^ 1) ^ 1 if (cs[0] == 1) else s.and_(p, q)
                if len(cs) == 1 and cs[0] == 1 and False: pass
                return sm, cy
        t = s.xor(a, c)
        sm = s.xor(t, b)
        cy = s.xor(s.and_(t, s.xor(b, c)), c)
        return sm, cy

    def add(s, x, y):
        c = F; out = []
        for i in range(32):
            sm, c = s.fa(x[i], y[i], c)
            out.append(sm)
        return out

    def addc(s, x, k): return s.add(x, const_word(k))

    def S0(s, x): return s.xorw(rotr(x, 2), rotr(x, 13), rotr(x, 22))
    def S1(s, x): return s.xorw(rotr(x, 6), rotr(x, 11), rotr(x, 25))
    def s0(s, x): return s.xorw(rotr(x, 7), rotr(x, 18), shr(x, 3))
    def s1(s, x): return s.xorw(rotr(x, 17), rotr(x, 19), shr(x, 10))

    def ch(s, e, f, g):
        return [s.xor(g[i], s.and_(e[i], s.xor(f[i], g[i]))) for i in range(32)]

    def maj(s, a, b, c):
        return [s.xor(b[i], s.and_(s.xor(a[i], b[i]), s.xor(b[i], c[i]))) for i in range(32)]


def inputs_words(x, nw, first=0):
    return [[x.input(first + 32 * k + 31 - j) for j in range(32)] for k in range(nw)]


def finish_spec(x, words):
    x.outputs = [w[31 - j] for w in words for j in range(32)]
    return x


def is_const(w): return all(v in (0, 1) for v in w)
def const_val(w): return sum(v << i for i, v in enumerate(w))
def nconst(w): return sum(1 for v in w if v in (0, 1))


def sumw(b, rows, consts=()):
    """Sum rows mod 2^32: all-constant rows and ints in consts folded first,
    variable rows added most-constant-bits first, const added to first var row."""
    k = sum(consts) & M
    var = []
    for r in rows:
        if is_const(r): k = (k + const_val(r)) & M
        else: var.append(r)
    if not var: return const_word(k)
    var.sort(key=lambda w: -nconst(w))
    acc = b.addc(var[0], k) if k else var[0]
    for r in var[1:]: acc = b.add(acc, r)
    return acc


def round_fn2(b, st, W, Kt):
    a, bb, c, d, e, f, g, h = st
    t1 = sumw(b, [h, b.S1(e), b.ch(e, f, g), W], (Kt,))
    t2 = sumw(b, [b.S0(a), b.maj(a, bb, c)])
    return [sumw(b, [t1, t2]), a, bb, c, sumw(b, [d, t1]), e, f, g]


def round_fn(b, st, W, Kt):
    a, bb, c, d, e, f, g, h = st
    t1 = b.add(b.add(b.add(h, b.S1(e)), b.ch(e, f, g)), b.addc(W, Kt))
    t2 = b.add(b.S0(a), b.maj(a, bb, c))
    return [b.add(t1, t2), a, bb, c, b.add(d, t1), e, f, g]


def build_round(t=0, round_f=round_fn):
    x = XAG(288); b = B(x)
    st = inputs_words(x, 9)
    out = round_f(b, st[:8], st[8], ref.K[t])
    return finish_spec(x, out)


def schedule_w(b, W):
    W = list(W)
    for t in range(16, 64):
        W.append(b.add(b.add(b.s1(W[t-2]), W[t-7]), b.add(b.s0(W[t-15]), W[t-16])))
    return W


def build_primary(round_f=round_fn, sched=schedule_w):
    x = XAG(512); b = B(x)
    Win = inputs_words(x, 16)
    W = sched(b, Win)
    st = [const_word(v) for v in ref.IV]
    for t in range(64):
        st = round_f(b, st, W[t], ref.K[t])
    out = [b.add(const_word(ref.IV[k]), st[k]) for k in range(8)]
    return finish_spec(x, out)


def report(x):
    y = cost.canonicalize(x)
    return y, cost.report_xag(y)


if __name__ == "__main__":
    import json
    what = sys.argv[1] if len(sys.argv) > 1 else "round"
    x = build_round(0) if what == "round" else build_primary()
    y, r = report(x)
    print(json.dumps(r))
