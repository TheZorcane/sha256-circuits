"""Faithful Python transliteration of emp-toolkit/emp-tool
emp-tool/circuits/crypto/sha256.h at merge commit 1623cb11 (PR #216, merged
2026-11-10 into prerelease/fixes), emitting a trusted XAG so the circuit can be
scored under the this work SPEC (NOT free, canonical AND+XOR) and simulated.

Bits are ('k', 0/1) public constants or ('w', lit) wires. Words are LSB-first
lists of 32 bits (emp convention: index i = bit i).

usage: python3 emp_sha256.py [var|fix] OUT.xag
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
from xag import XAG  # noqa: E402
import sha256_ref as ref  # noqa: E402


class Ctx:
    def __init__(self, n):
        self.x = XAG(n)

    def xor_gate(self, a, b):
        return self.x.XOR(a, b)

    def and_gate(self, a, b):
        return self.x.AND(a, b)

    def not_gate(self, a):
        return a ^ 1


def kbit(v):
    return ("k", 1 if v else 0)


def wbit(l):
    return ("w", l)


def bnot(c, a):
    return kbit(a[1] == 0) if a[0] == "k" else wbit(c.not_gate(a[1]))


def bxor(c, a, b):
    if a[0] == "k" and b[0] == "k":
        return kbit(a[1] != b[1])
    if a[0] == "k":
        return bnot(c, b) if a[1] else b
    if b[0] == "k":
        return bnot(c, a) if b[1] else a
    return wbit(c.xor_gate(a[1], b[1]))


def band(c, a, b):
    if (a[0] == "k" and a[1] == 0) or (b[0] == "k" and b[1] == 0):
        return kbit(0)
    if a[0] == "k":
        return b
    if b[0] == "k":
        return a
    return wbit(c.and_gate(a[1], b[1]))


def kconst(v):
    return [kbit((v >> i) & 1) for i in range(32)]


def all_known(x):
    return all(b[0] == "k" for b in x)


def any_wire(x):
    return not all_known(x)


def ksig(c, x, r1, r2, r3, shift3):
    r = []
    for i in range(32):
        t3 = kbit(0) if (shift3 and i + r3 >= 32) else x[(i + r3) % 32]
        r.append(bxor(c, bxor(c, x[(i + r1) % 32], x[(i + r2) % 32]), t3))
    return r


def kch(c, e, f, g):
    r = []
    for i in range(32):
        if e[i][0] == "k":
            r.append(f[i] if e[i][1] else g[i])
            continue
        r.append(bxor(c, g[i], band(c, e[i], bxor(c, f[i], g[i]))))
    return r


def kmaj(c, a, b, cc):
    r = []
    for i in range(32):
        x, y, z = a[i], b[i], cc[i]
        if x[0] == "k" and y[0] == "k":
            r.append(x if x[1] == y[1] else z); continue
        if x[0] == "k" and z[0] == "k":
            r.append(x if x[1] == z[1] else y); continue
        if y[0] == "k" and z[0] == "k":
            r.append(y if y[1] == z[1] else x); continue
        r.append(bxor(c, x, band(c, bxor(c, x, y), bxor(c, x, z))))
    return r


def ksum(c, terms, k):
    k &= 0xFFFFFFFF
    col = [[] for _ in range(33)]
    for t in terms:
        for i in range(32):
            b = t[i]
            if b[0] == "k":
                k = (k + (b[1] << i)) & 0xFFFFFFFF
            else:
                col[i].append(b[1])
    r = [None] * 32
    for i in range(32):
        v = col[i]
        kb = (k >> i) & 1
        if i == 31:
            s = kbit(kb)
            for w in v:
                s = bxor(c, s, wbit(w))
            r[i] = s
            break
        p = 0
        while len(v) - p >= 3:
            x, y, z = v[p], v[p + 1], v[p + 2]
            p += 3
            xy = c.xor_gate(x, y)
            v.append(c.xor_gate(xy, z))
            col[i + 1].append(c.xor_gate(x, c.and_gate(xy, c.xor_gate(x, z))))
        left = len(v) - p
        if left == 0:
            r[i] = kbit(kb)
        elif left == 1:
            x = v[p]
            if not kb:
                r[i] = wbit(x)
            else:
                r[i] = wbit(c.not_gate(x))
                col[i + 1].append(x)
        else:
            x, y = v[p], v[p + 1]
            xy, a = c.xor_gate(x, y), c.and_gate(x, y)
            if not kb:
                r[i] = wbit(xy)
                col[i + 1].append(a)
            else:
                r[i] = wbit(c.not_gate(xy))
                col[i + 1].append(c.xor_gate(xy, a))
    return r


def compress_k(c, state, block):
    W = list(block) + [None] * 48
    for t in range(16, 64):
        s1 = ksig(c, W[t - 2], 17, 19, 10, True)
        s0 = ksig(c, W[t - 15], 7, 18, 3, True)
        W[t] = ksum(c, [s1, W[t - 7], s0, W[t - 16]], 0)
    a, b, cc, d, e, f, g, h = state
    for t in range(64):
        S1 = ksig(c, e, 6, 11, 25, False)
        CH = kch(c, e, f, g)
        S0 = ksig(c, a, 2, 13, 22, False)
        MJ = kmaj(c, a, b, cc)
        t1w = any_wire(h) + any_wire(S1) + any_wire(CH) + any_wire(W[t])
        if t1w <= 1 and all_known(d) and all_known(S0) and all_known(MJ):
            e_new = ksum(c, [h, W[t], S1, CH, d], ref.K[t])
            a_new = ksum(c, [h, W[t], S1, CH, S0, MJ], ref.K[t])
        else:
            T1 = ksum(c, [h, W[t], S1, CH], ref.K[t])
            e_new = ksum(c, [d, T1], 0)
            a_new = ksum(c, [S0, MJ, T1], 0)
        h, g, f, e = g, f, e, e_new
        d, cc, b, a = cc, b, a, a_new
    v = [a, b, cc, d, e, f, g, h]
    return [ksum(c, [state[i], v[i]], 0) for i in range(8)]


def build(mode):
    """SPEC order: SPEC bit j of word w = bit 31-j (MSB first)."""
    n = 768 if mode == "var" else 512
    c = Ctx(n)
    x = c.x

    def word_in(base_bit):
        return [wbit(x.input(base_bit + 31 - i)) for i in range(32)]  # LSB-first

    if mode == "var":
        state = [word_in(32 * w) for w in range(8)]
        block = [word_in(256 + 32 * w) for w in range(16)]
    else:
        state = [kconst(v) for v in ref.IV]
        block = [word_in(32 * w) for w in range(16)]
    out = compress_k(c, state, block)
    outs = []
    for w in range(8):
        for j in range(32):
            b = out[w][31 - j]
            outs.append(b[1] if b[0] == "w" else b[1])  # constant lit 0/1
    x.outputs = outs
    x.validate()
    return x


if __name__ == "__main__":
    from xag import write_xag
    x = build(sys.argv[1])
    open(sys.argv[2], "w").write(write_xag(x))
    print(sys.argv[1], x.counts())
