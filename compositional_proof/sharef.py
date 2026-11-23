"""Independent word-level -> XAG reference generators for SHA-256 (FIPS 180-4).
A word is a list of 32 literals, LSB first (index j = bit j).  Constants are lits 0/1.
Two styles: style 'A' (textbook, ripple adders) and style 'B' (alternative: Kogge-Stone
adders, sigma via GF(2) matrices derived from trusted sha256_ref, ch/maj in mux form)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
import sha256_ref as S
from xag import XAG


def cword(v): return [(v >> j) & 1 for j in range(32)]
def rotr(w, n): return [w[(j + n) % 32] for j in range(32)]
def shr(w, n): return [w[j + n] if j + n < 32 else 0 for j in range(32)]
def bxor(x, a, b): return [x.XOR(p, q) for p, q in zip(a, b)]


# ---------------- style A (textbook) ----------------
def add_A(x, a, b):
    c, out = 0, []
    for j in range(32):
        p = x.XOR(a[j], b[j])
        out.append(x.XOR(p, c))
        if j < 31:
            c = x.XOR(x.AND(a[j], b[j]), x.AND(c, p))
    return out

def x3(x, a, b, c): return bxor(x, bxor(x, a, b), c)
def bsig0_A(x, w): return x3(x, rotr(w, 2), rotr(w, 13), rotr(w, 22))
def bsig1_A(x, w): return x3(x, rotr(w, 6), rotr(w, 11), rotr(w, 25))
def ssig0_A(x, w): return x3(x, rotr(w, 7), rotr(w, 18), shr(w, 3))
def ssig1_A(x, w): return x3(x, rotr(w, 17), rotr(w, 19), shr(w, 10))
def ch_A(x, e, f, g):  # (e&f) ^ (~e&g)
    return [x.XOR(x.AND(e[j], f[j]), x.AND(e[j] ^ 1, g[j])) for j in range(32)]
def maj_A(x, a, b, c):  # (a&b)^(a&c)^(b&c)
    return [x.XOR(x.XOR(x.AND(a[j], b[j]), x.AND(a[j], c[j])), x.AND(b[j], c[j])) for j in range(32)]


# ---------------- style B (alternative construction) ----------------
def add_B(x, a, b):  # Kogge-Stone
    n = 32
    p = [x.XOR(a[i], b[i]) for i in range(n)]
    g = [x.AND(a[i], b[i]) for i in range(n)]
    G, P = g[:], p[:]
    d = 1
    while d < n:
        G2, P2 = G[:], P[:]
        for i in range(d, n):
            G2[i] = x.XOR(G[i], x.AND(P[i], G[i - d]))   # disjoint -> OR==XOR
            P2[i] = x.AND(P[i], P[i - d])
        G, P = G2, P2
        d *= 2
    return [p[0]] + [x.XOR(p[i], G[i - 1]) for i in range(1, n)]

def _matrix(f):
    """GF(2) matrix of word function f taken from trusted sha256_ref (unit-vector probing)."""
    rows = [[] for _ in range(32)]
    for j in range(32):
        y = f(1 << j)
        for i in range(32):
            if (y >> i) & 1:
                rows[i].append(j)
    import random
    r = random.Random(5)
    for _ in range(50):
        u, v = r.getrandbits(32), r.getrandbits(32)
        assert f(u ^ v) == f(u) ^ f(v)
    return rows
MAT = {"bsig0": _matrix(S._bsig0), "bsig1": _matrix(S._bsig1),
       "ssig0": _matrix(S._ssig0), "ssig1": _matrix(S._ssig1)}
def _lin(name):
    def f(x, w):
        out = []
        for row in MAT[name]:
            l = w[row[0]]
            for j in row[1:]:
                l = x.XOR(l, w[j])
            out.append(l)
        return out
    return f
bsig0_B, bsig1_B, ssig0_B, ssig1_B = (_lin(n) for n in ("bsig0", "bsig1", "ssig0", "ssig1"))
def ch_B(x, e, f, g):  # g ^ (e & (f^g))
    return [x.XOR(g[j], x.AND(e[j], x.XOR(f[j], g[j]))) for j in range(32)]
def maj_B(x, a, b, c):  # a ^ ((a^b)&(a^c))
    return [x.XOR(a[j], x.AND(x.XOR(a[j], b[j]), x.XOR(a[j], c[j]))) for j in range(32)]

STY = {"A": dict(add=add_A, bsig0=bsig0_A, bsig1=bsig1_A, ssig0=ssig0_A, ssig1=ssig1_A, ch=ch_A, maj=maj_A),
       "B": dict(add=add_B, bsig0=bsig0_B, bsig1=bsig1_B, ssig0=ssig0_B, ssig1=ssig1_B, ch=ch_B, maj=maj_B)}


# ---------------- components ----------------
def sched_step(x, w2, w7, w15, w16, st="A"):
    """W_t = ssig1(W_{t-2}) + W_{t-7} + ssig0(W_{t-15}) + W_{t-16}"""
    s = STY[st]
    r = s["add"](x, s["ssig1"](x, w2), w7)
    r = s["add"](x, r, s["ssig0"](x, w15))
    return s["add"](x, r, w16)

def round_step(x, state, W, K, st="A"):
    """state=(a,b,c,d,e,f,g,h); returns (a',e') (the two new words; the rest are shifts)."""
    a, b, c, d, e, f, g, h = state
    s = STY[st]
    t1 = s["add"](x, h, s["bsig1"](x, e))
    t1 = s["add"](x, t1, s["ch"](x, e, f, g))
    t1 = s["add"](x, t1, cword(K))
    t1 = s["add"](x, t1, W)
    t2 = s["add"](x, s["bsig0"](x, a), s["maj"](x, a, b, c))
    return s["add"](x, t1, t2), s["add"](x, d, t1)

def final_add(x, iv, w, st="A"):
    return STY[st]["add"](x, cword(iv), w)


# ---------------- word-level wiring (single source of truth) ----------------
# symbolic word names: ("W",t) t=0..63 ; ("a",t),("e",t) t=-4..63 (t<0: IV constants)
def sched_inputs(t): return [("W", t - 2), ("W", t - 7), ("W", t - 15), ("W", t - 16)]
def round_inputs(t):
    return [("a", t - 1), ("a", t - 2), ("a", t - 3), ("a", t - 4),
            ("e", t - 1), ("e", t - 2), ("e", t - 3), ("e", t - 4), ("W", t)]
def final_inputs(k):  # H_out[k] = IV[k] + (a_{63-k} for k<4 ; e_{63-(k-4)} for k>=4)
    return [("a", 63 - k)] if k < 4 else [("e", 63 - (k - 4))]
def iv_word(name):
    kind, t = name
    if t >= 0:
        return None
    i = -t - 1  # a_{-1}=IV0,a_{-2}=IV1,a_{-3}=IV2,a_{-4}=IV3 ; e_{-1}=IV4,...
    return S.IV[i] if kind == "a" else S.IV[4 + i]

# column-compression generator modification (copy only, env VERIFY_TAIL=1): sched 62/63, round 62/63 and final 0..7 replaced by
# one "tail" lemma whose cone starts at the round-62 state words and the schedule sources of W62/W63.
import os as _os
TAIL = _os.environ.get("VERIFY_TAIL") == "1"
TAIL_IN = [("a", 61), ("a", 60), ("a", 59), ("a", 58), ("e", 61), ("e", 60), ("e", 59), ("e", 58),
           ("W", 60), ("W", 55), ("W", 47), ("W", 46), ("W", 61), ("W", 56), ("W", 48)]
SPLIT = _os.environ.get("VERIFY_SPLIT") == "1"     # carry-root generator modification (copy only)
R62_IN = [("a", 61), ("a", 60), ("a", 59), ("a", 58), ("e", 61), ("e", 60), ("e", 59), ("e", 58),
          ("W", 60), ("W", 55), ("W", 47), ("W", 46)]
T2_IN = [("a", 62), ("a", 61), ("a", 60), ("a", 59), ("e", 62), ("e", 61), ("e", 60), ("e", 59),
         ("W", 61), ("W", 56), ("W", 48), ("W", 47)]
def r62f_apply(x, inw, st="A"):
    w62 = sched_step(x, inw[("W", 60)], inw[("W", 55)], inw[("W", 47)], inw[("W", 46)], st=st)
    a2, e2 = round_step(x, [inw[n] for n in R62_IN[:8]], w62, S.K[62], st=st)
    return [a2, e2]
def tail2_apply(x, inw, st="A"):
    w63 = sched_step(x, inw[("W", 61)], inw[("W", 56)], inw[("W", 48)], inw[("W", 47)], st=st)
    st8 = [inw[n] for n in T2_IN[:8]]
    a2, e2 = round_step(x, st8, w63, S.K[63], st=st)
    st8 = [a2] + st8[:3] + [e2] + st8[4:7]
    return [final_add(x, S.IV[k], st8[k], st=st) for k in range(8)]
def lemma_list():
    if TAIL and SPLIT:
        return ([("sched", t) for t in range(16, 62)] + [("round", t) for t in range(62)] + [("r62f", 0), ("tail2", 0)])
    if TAIL:
        return ([("sched", t) for t in range(16, 62)] + [("round", t) for t in range(62)] + [("tail", 0)])
    return ([("sched", t) for t in range(16, 64)] + [("round", t) for t in range(64)]
            + [("final", k) for k in range(8)])
def tail_apply(x, inw, st="A"):
    W = dict((n[1], inw[n]) for n in TAIL_IN if n[0] == "W")
    W[62] = sched_step(x, W[60], W[55], W[47], W[46], st=st)
    W[63] = sched_step(x, W[61], W[56], W[48], W[47], st=st)
    st8 = [inw[n] for n in TAIL_IN[:8]]
    for t in (62, 63):
        a2, e2 = round_step(x, st8, W[t], S.K[t], st=st)
        st8 = [a2] + st8[:3] + [e2] + st8[4:7]
    return [final_add(x, S.IV[k], st8[k], st=st) for k in range(8)]
def lemma_in_out(l):
    kind, i = l
    if kind == "tail": return list(TAIL_IN), [("out", k) for k in range(8)]
    if kind == "r62f": return list(R62_IN), [("a", 62), ("e", 62)]
    if kind == "tail2": return list(T2_IN), [("out", k) for k in range(8)]
    if kind == "sched": return sched_inputs(i), [("W", i)]
    if kind == "round": return round_inputs(i), [("a", i), ("e", i)]
    return final_inputs(i), [("out", i)]
def lemma_apply(x, l, inw, st="A"):
    """inw: dict name->word (lits). returns list of output words."""
    kind, i = l
    if kind == "tail":
        return tail_apply(x, inw, st)
    if kind == "r62f":
        return r62f_apply(x, inw, st)
    if kind == "tail2":
        return tail2_apply(x, inw, st)
    if kind == "sched":
        return [sched_step(x, *[inw[k] for k in sched_inputs(i)], st=st)]
    if kind == "round":
        n = round_inputs(i)
        return list(round_step(x, [inw[k] for k in n[:8]], inw[n[8]], S.K[i], st=st))
    return [final_add(x, S.IV[i], inw[final_inputs(i)[0]], st=st)]

def full_primary(st="A"):
    """Full reference PRIMARY (512 -> 256) built by chaining lemma_apply over the wiring."""
    x = XAG(512)
    vals = {}
    for t in range(16):
        vals[("W", t)] = [x.input(32 * t + 31 - j) for j in range(32)]
    for nm in [(k, t) for k in "ae" for t in range(-4, 0)]:
        vals[nm] = cword(iv_word(nm))
    outs = [None] * 8
    for l in lemma_list():
        ins, outn = lemma_in_out(l)
        res = lemma_apply(x, l, vals, st)
        if l[0] == "final":
            outs[l[1]] = res[0]
        else:
            for nm, w in zip(outn, res):
                vals[nm] = w
    x.outputs = [outs[k][31 - j] for k in range(8) for j in range(32)]
    return x
