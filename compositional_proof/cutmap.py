"""Cut map: word slot (name, j) -> list of candidate literals whose value == that bit of the word
for ALL inputs (to be proven by the lemma that produces the word).
 - W_t (t<16): the primary input literal (by definition of the SPEC bit order), only alternative.
 - ("out",k): the actual primary output literal, only alternative.
 - others: all nodes matching the reference-bit simulation signature (cuts.py)."""
import common, json
from xag import read_xag


def load(path="cand.xag"):
    x = read_xag(open(path).read())
    m = json.load(open("work/cutmap.json"))
    return x, m


def build_cutmap(x, m):
    """returns dict (name,j) -> list of lits (value(lit) == word bit)"""
    import cuts, random
    import sha256_ref as S
    # re-derive all matches (not only the first) from simulation
    N = cuts.N
    r = random.Random(20261118)
    vs = [[r.getrandbits(32) for _ in range(16)] for _ in range(N)]
    vs[0] = [0] * 16; vs[1] = [S.MASK] * 16
    import sim
    lns_ = sim.vectors_to_lanes([S.words_to_int(w) for w in vs], 512)
    val = cuts.node_values(x, lns_, N)
    mask = (1 << N) - 1
    sig = {}
    for v in range(1, len(val)):
        sig.setdefault(val[v], []).append(2 * v)             # positive literal has this sig
        sig.setdefault(val[v] ^ mask, []).append(2 * v + 1)  # complemented literal
    per = [cuts.ref_words(w) for w in vs]
    cm = {}
    for t in range(64):
        for kind, idx in (("W", 0), ("a", 1), ("e", 2)):
            for j in range(32):
                s = 0
                for k, p in enumerate(per):
                    s |= ((p[idx][t] if idx == 0 else p[idx][t]) >> j & 1) << k
                if kind == "W" and t < 16:
                    cm[(("W", t), j)] = [2 * (32 * t + 31 - j + 1)]
                else:
                    cm[((kind, t), j)] = sorted(sig.get(s, []))
    for k in range(8):
        for j in range(32):
            cm[(("out", k), j)] = [x.outputs[32 * k + 31 - j]]
    # primary-only: keep the single most-used match (ties: lowest literal); other matches are only
    # signature coincidences or unused duplicates and are treated as ordinary interior nodes.
    fan = {}
    for op, a, b in x.nodes:
        for l in (a, b):
            fan[l >> 1] = fan.get(l >> 1, 0) + 1
    ALL = {}
    for k, v in cm.items():
        if len(v) > 1:
            ALL[k] = list(v)
            cm[k] = [sorted(v, key=lambda L: (-fan.get(L >> 1, 0), L))[0]]
    global DROPPED
    DROPPED = ALL
    return cm
