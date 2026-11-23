"""port generator SHA-1 compression generator (fixed IV, 512 -> 160, SPEC in sha1ref.py) using the column-compression generator column-DP multi-operand
adder (madd) : one 4-row + constant sum per round, Ch-pair trick in rounds 0..19, feed-forward IV0 fused into round 79.
usage: gen_sha1.py OUT.xag [variant]   tokens: pair (Ch pair), fuse (IV0 fused into last round)"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gen3
from gen3 import XAG, X, A, cost, PC, write_xag
import madd
import sha1ref as R


def build(variant="pair+fuse"):
    ws = 32; M = R.M
    x = XAG(512)
    W = [[x.input(32 * k + 31 - j) for j in range(32)] for k in range(16)]
    def S(rows, K=0, pairs=()): return madd.build_sum(x, rows, K & M, n=ws, plan_cache=PC, pairs=pairs)
    def rotl(w, n): return [w[(i - n) % 32] for i in range(32)]
    def xorw(*w_):
        r = w_[0]
        for w in w_[1:]: r = [X(x, p, q) for p, q in zip(r, w)]
        return r
    for t in range(16, 80):
        W.append(rotl(xorw(W[t-3], W[t-8], W[t-14], W[t-16]), 1))
    cw = lambda v: [(v >> i) & 1 for i in range(32)]
    a, b, c, d, e = [cw(v) for v in R.IV]
    for t in range(80):
        Kt = R.K[t // 20]
        if t < 20:
            # Ch(b,c,d) = d ^ (b & (c^d));  e ^ Ch = (e^d) ^ (b & (c^d)); e^d is round t-1's c^d (strash-shared)
            if "pair" in variant:
                pr = [X(x, X(x, d[i], e[i]), A(x, b[i], X(x, c[i], d[i]))) for i in range(32)]
                rows, pairs = [rotl(a, 5), W[t]], [(e, pr)]
            else:
                f = [X(x, d[i], A(x, b[i], X(x, c[i], d[i]))) for i in range(32)]
                rows, pairs = [rotl(a, 5), f, e, W[t]], []
        elif t < 40 or t >= 60:
            rows, pairs = [rotl(a, 5), xorw(b, c, d), e, W[t]], []
        else:
            mj = [X(x, c[i], A(x, X(x, b[i], c[i]), X(x, c[i], d[i]))) for i in range(32)]
            rows, pairs = [rotl(a, 5), mj, e, W[t]], []
        K2 = Kt + (R.IV[0] if (t == 79 and "fuse" in variant) else 0)
        T = S(rows, K2, pairs)
        a, b, c, d, e = T, a, rotl(b, 30), c, d
    st = [a, b, c, d, e]
    out = [a if "fuse" in variant else S([a], R.IV[0])] + [S([st[k]], R.IV[k]) for k in range(1, 5)]
    x.outputs = [w[31 - j] for w in out for j in range(32)]
    return x


if __name__ == "__main__":
    variant = sys.argv[2] if len(sys.argv) > 2 else "pair+fuse"
    x = build(variant)
    y = cost.canonicalize(x)
    print(json.dumps(cost.report_xag(y)))
    open(sys.argv[1], "w").write(write_xag(y))
