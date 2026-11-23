import random, sys
from xag2 import parse
def run(path, N, vecs):
    n, nodes, outs = parse(path); M = (1 << N) - 1
    v = [0] + vecs
    for k, a, b in nodes:
        x = v[a >> 1] ^ (M if a & 1 else 0); y = v[b >> 1] ^ (M if b & 1 else 0)
        v.append(x & y if k == 'A' else x ^ y)
    return [v[o >> 1] ^ (M if o & 1 else 0) for o in outs]
if __name__ == "__main__":
    ref, new = sys.argv[1], sys.argv[2]; N = int(sys.argv[3])
    rnd = random.Random(12345)
    n = parse(ref)[0]
    vecs = [rnd.getrandbits(N) for _ in range(n)]
    r = run(ref, N, vecs); w = run(new, N, vecs)
    bad = [i for i in range(len(r)) if r[i] != w[i]]
    print("vectors", N, "outputs", len(r), "mismatching outputs", len(bad), bad[:5])
    sys.exit(1 if bad else 0)
