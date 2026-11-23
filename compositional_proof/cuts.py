"""Step 2: independent cut identification by simulation signatures against the verifier's own reference
semantics (word values computed with verifier/sha256_ref primitives), then a cut map
  word name -> list of 32 candidate literals (LSB first), const lits 0/1 allowed.
Writes work/cutmap.json."""
import common, json, random, sys
import sha256_ref as S
from xag import read_xag, AND
import cost, sim
from sharef import iv_word

N = 16384


def node_values(x, lns_, n):
    mask = (1 << n) - 1
    val = [0] + [v & mask for v in lns_]
    for op, a, b in x.nodes:
        va = val[a >> 1] ^ (mask if a & 1 else 0)
        vb = val[b >> 1] ^ (mask if b & 1 else 0)
        val.append(va & vb if op == AND else va ^ vb)
    return val


def ref_words(w16):
    """all intermediate words for one block, using trusted sha256_ref primitives"""
    W = list(w16)
    for t in range(16, 64):
        W.append((S._ssig1(W[t - 2]) + W[t - 7] + S._ssig0(W[t - 15]) + W[t - 16]) & S.MASK)
    a, b, c, d, e, f, g, h = S.IV
    A, E = {}, {}
    for t in range(64):
        t1 = (h + S._bsig1(e) + S._ch(e, f, g) + S.K[t] + W[t]) & S.MASK
        t2 = (S._bsig0(a) + S._maj(a, b, c)) & S.MASK
        h, g, f, e, d, c, b, a = g, f, e, (d + t1) & S.MASK, c, b, a, (t1 + t2) & S.MASK
        A[t], E[t] = a, e
    out = [(x + y) & S.MASK for x, y in zip(S.IV, (a, b, c, d, e, f, g, h))]
    return W, A, E, out


def main(path="cand.xag"):
    x = read_xag(open(path).read())
    r = random.Random(20261118)
    vs = [[r.getrandbits(32) for _ in range(16)] for _ in range(N)]
    vs[0] = [0] * 16; vs[1] = [S.MASK] * 16
    lns_ = sim.vectors_to_lanes([S.words_to_int(w) for w in vs], 512)
    val = node_values(x, lns_, N)
    mask = (1 << N) - 1
    sig = {}
    for v in range(1, len(val)):
        s = val[v]
        sig.setdefault(s, []).append(v)
        sig.setdefault(s ^ mask, []).append(-v)   # negative: complemented
    words = {}
    per = [ref_words(w) for w in vs]
    for t in range(64):
        words[("W", t)] = [p[0][t] for p in per]
        words[("a", t)] = [p[1][t] for p in per]
        words[("e", t)] = [p[2][t] for p in per]
    for k in range(8):
        words[("out", k)] = [p[3][k] for p in per]
    cut, multi, missing = {}, {}, []
    for name, vals in words.items():
        lits = []
        for j in range(32):
            s = 0
            for k, w in enumerate(vals):
                s |= ((w >> j) & 1) << k
            c = sig.get(s)
            if not c:
                missing.append((name, j)); lits.append(None); continue
            # prefer positive polarity, earliest node
            c = sorted(c, key=lambda z: (z < 0, abs(z)))
            if len(c) > 1:
                multi[f"{name}:{j}"] = c
            z = c[0]
            lits.append(2 * abs(z) + (1 if z < 0 else 0))   # literal such that value == word bit
        cut[name] = lits
    print("missing", missing[:10], len(missing)); print("multi", len(multi), list(multi.items())[:5])
    json.dump({"cut": {f"{k[0]}{k[1]}": v for k, v in cut.items()}, "multi": multi, "missing": missing,
               "cand_hash": cost.xag_hash(cost.canonicalize(x)), "n_vectors": N},
              open("work/cutmap.json", "w"))
    return x, cut, multi


if __name__ == "__main__":
    main()
