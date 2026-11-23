"""Plain 'textbook' circuits as calibration baselines.

Chained 2-operand ripple-carry adders (1-AND full adder: s=a^b^c, c'=c^((a^c)&(b^c))),
Ch = g^(e&(f^g)), Maj = b^((a^b)&(b^c)), parity f = b^c^d, sigma/Sigma as XORs.
Constants (IV, K) are wired as constants; trusted canonicalize folds them.
No MDFA, no carry-save, no cross-round sharing.  SPEC-style bit order
(words in order, MSB first).

usage: python3 textbook.py OUT.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
sys.path.insert(0, HERE)
from xag import XAG  # noqa: E402
import cost  # noqa: E402
import sim  # noqa: E402
import sha256_ref as ref256  # noqa: E402
import refs  # noqa: E402
from rescore import verify_generic, not_tally  # noqa: E402


class B:
    def __init__(self, n, w):
        self.x, self.w = XAG(n), w

    def inword(self, base):  # LSB-first list of lits; SPEC bit j of word = bit w-1-j
        return [self.x.input(base + self.w - 1 - i) for i in range(self.w)]

    def const(self, v):
        return [(v >> i) & 1 for i in range(self.w)]  # lit 0 / 1

    def xor(self, a, b):
        return [self.x.XOR(p, q) for p, q in zip(a, b)]

    def rotr(self, a, r):
        return [a[(i + r) % self.w] for i in range(self.w)]

    def shr(self, a, r):
        return [a[i + r] if i + r < self.w else 0 for i in range(self.w)]

    def add(self, a, b):
        X = self.x
        out, c = [], None
        for i in range(self.w):
            if c is None:
                out.append(X.XOR(a[i], b[i]))
                c = X.AND(a[i], b[i])
                continue
            ac, bc = X.XOR(a[i], c), X.XOR(b[i], c)
            out.append(X.XOR(ac, b[i]))
            if i < self.w - 1:
                c = X.XOR(c, X.AND(ac, bc))
        return out

    def sum(self, *ws):
        acc = ws[0]
        for t in ws[1:]:
            acc = self.add(acc, t)
        return acc

    def ch(self, e, f, g):
        X = self.x
        return [X.XOR(gi, X.AND(ei, X.XOR(fi, gi))) for ei, fi, gi in zip(e, f, g)]

    def maj(self, a, b, c):
        X = self.x
        return [X.XOR(bi, X.AND(X.XOR(ai, bi), X.XOR(bi, ci))) for ai, bi, ci in zip(a, b, c)]

    def finish(self, words):
        self.x.outputs = [w[self.w - 1 - j] for w in words for j in range(self.w)]
        return self.x


def sha2(w, fixed):
    if w == 32:
        K, IV = ref256.K, ref256.IV
        s0, s1, S0, S1, R = (7, 18, 3), (17, 19, 10), (2, 13, 22), (6, 11, 25), 64
    else:
        K, IV = refs.SHA512_K, refs.SHA512_IV
        s0, s1, S0, S1, R = (1, 8, 7), (19, 61, 6), (28, 34, 39), (14, 18, 41), 80
    b = B(16 * w if fixed else 24 * w, w)
    if fixed:
        cv = [b.const(v) for v in IV]
        W = [b.inword(w * i) for i in range(16)]
    else:
        cv = [b.inword(w * i) for i in range(8)]
        W = [b.inword(8 * w + w * i) for i in range(16)]
    sig = lambda x, r: b.xor(b.xor(b.rotr(x, r[0]), b.rotr(x, r[1])), b.shr(x, r[2]))  # noqa: E731
    Sig = lambda x, r: b.xor(b.xor(b.rotr(x, r[0]), b.rotr(x, r[1])), b.rotr(x, r[2]))  # noqa: E731
    for t in range(16, R):
        W.append(b.sum(sig(W[t - 2], s1), W[t - 7], sig(W[t - 15], s0), W[t - 16]))
    a, bb, c, d, e, f, g, h = cv
    for t in range(R):
        T1 = b.sum(h, Sig(e, S1), b.ch(e, f, g), b.const(K[t]), W[t])
        T2 = b.add(Sig(a, S0), b.maj(a, bb, c))
        h, g, f, e, d, c, bb, a = g, f, e, b.add(d, T1), c, bb, a, b.add(T1, T2)
    return b.finish([b.add(x, y) for x, y in zip(cv, (a, bb, c, d, e, f, g, h))])


def sha1():
    b = B(512, 32)
    cv = [b.const(v) for v in refs.SHA1_IV]
    W = [b.inword(32 * i) for i in range(16)]
    for t in range(16, 80):
        x = b.xor(b.xor(W[t - 3], W[t - 8]), b.xor(W[t - 14], W[t - 16]))
        W.append(b.rotr(x, 31))
    a, bb, c, d, e = cv
    for t in range(80):
        if t < 20:
            f, k = b.ch(bb, c, d), 0x5A827999
        elif t < 40:
            f, k = b.xor(b.xor(bb, c), d), 0x6ED9EBA1
        elif t < 60:
            f, k = b.maj(bb, c, d), 0x8F1BBCDC
        else:
            f, k = b.xor(b.xor(bb, c), d), 0xCA62C1D6
        tmp = b.sum(b.rotr(a, 27), f, e, b.const(k), W[t])
        a, bb, c, d, e = tmp, a, b.rotr(bb, 2), c, d
    return b.finish([b.add(x, y) for x, y in zip(cv, (a, bb, c, d, e))])


if __name__ == "__main__":
    res = {}
    for name, mk, kind in [("textbook_sha256_fix", lambda: sha2(32, True), "primary"),
                           ("textbook_sha256_var", lambda: sha2(32, False), "secondary"),
                           ("textbook_sha1_fix", sha1, "sha1_fix"),
                           ("textbook_sha512_var", lambda: sha2(64, False), "sha512_var"),
                           ("textbook_sha512_fix", lambda: sha2(64, True), "sha512_fix")]:
        x = cost.canonicalize(mk())
        r = cost.report_xag(x)
        if kind in ("primary", "secondary"):
            v = sim.verify(x, kind, list(range(x.n_inputs)), [(i, 0) for i in range(256)], n=5000, seed=7)
            v = {k: v[k] for k in ("status", "vectors_tested", "mismatches")}
        else:
            v = verify_generic(x, kind, n=4000, seed=7)
        res[name] = {"canonical": {k: r[k] for k in ("AND", "XOR", "total", "and_depth")},
                     "sha256": r["sha256"], "not_tally": not_tally(x), "sim": v}
        print(name, json.dumps(res[name]), flush=True)
    json.dump(res, open(sys.argv[1], "w"), indent=1)
    bad = [k for k, v in res.items() if v["sim"].get("status") != "pass"]
    if bad:
        print("SIMULATION FAILED:", ", ".join(bad))
        sys.exit(1)
