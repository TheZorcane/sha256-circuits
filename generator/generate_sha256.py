"""SHA-256 generator for the 96,273 and 96,324 lines: carry-root choice with difference-word chains (F-chains),
self-F and self-delta by-products, combined with the AND-weighted planner objective and the joint-sum planner.
It extends the constant-redistribution generator gen4.py (cshift+neg+final+wfuse).
A schedule word W_j is a plain operand of up to three sums: T1_j, W_{j+7}, W_{j+16}.  With the difference word
F(W_j)[i] = W_j[i-1] ^ W_j[i] built once, each of those sums can run an F-chain along W_j (madd3.py) and save
one gate per column.
usage: generate_sha256.py OUT.xag [J]     J = 'all' | 'cover' | 'none' | comma list of word indices | 'auto'
env GEN_T1F=0 disables F-chains in T1 (schedule only); see README.md for the other switches."""
import sys, os, json, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen3, gen4
from gen3 import XAG, FB, g2, sref, cost, chpair, write_xag
import madd, madd3 as madd2

M32 = (1 << 32) - 1
import pickle, hashlib
_PCF = os.path.join(os.environ.get("GEN_CACHE", os.path.join(HERE, "cache")), "plancache_%s.pkl" % hashlib.sha256(
    (repr(sorted(madd2.ENABLED)) + repr(madd2.SELF) + repr(madd2.LAM) + repr(sorted(madd.OPS)) + open(os.path.join(HERE, "madd3.py")).read()).encode()).hexdigest()[:12])
try: PC2 = pickle.load(open(_PCF, "rb"))
except Exception: PC2 = {}
_PC0 = len(PC2)


def save_cache():
    if len(PC2) != _PC0:
        os.makedirs(os.path.dirname(_PCF), exist_ok=True)
        pickle.dump(PC2, open(_PCF, "wb"))
T1F = os.environ.get("GEN_T1F", "1") == "1"


def negate(rows, K, pairs, frows, dwire=None):
    n = len(rows) + 2 * len(pairs) + len(frows) + (1 if dwire else 0)
    r2 = [[l ^ 1 for l in r] for r in rows]
    p2 = [([u ^ 1 for u in U], Wd) for (U, Wd) in pairs]
    f2 = [[l ^ 1 for l in r] for r in frows]
    d2 = ([l ^ 1 for l in dwire[0]], dwire[1]) if dwire else None      # ~Maj = ~root ^ M
    return r2, (n - 1 - K) & M32, p2, f2, d2


D32 = [1] * 32


def E(rows, K, pairs=(), neg=False, frows=(), wantF=False, d=False):
    if neg: rows, K, pairs, frows, _ = negate(rows, K, pairs, frows, None if not d else ([2] * 32, [2] * 32))
    return madd2.estimate(rows, K & M32, pairs=pairs, frows=frows, plan_cache=PC2, wantF=wantF, dcols=D32 if d else None)


def S(x, rows, K, pairs=(), neg=False, frows=(), wantF=False, dwire=None):
    if neg:
        r2, K2, p2, f2, d2 = negate(rows, K, pairs, frows, dwire)
        return [l ^ 1 for l in madd2.build_sum(x, r2, K2, plan_cache=PC2, pairs=p2, frows=f2, XF=gen3.X, AF=gen3.A, wantF=wantF, dwire=d2)]
    return madd2.build_sum(x, rows, K & M32, plan_cache=PC2, pairs=pairs, frows=frows, XF=gen3.X, AF=gen3.A, wantF=wantF, dwire=dwire)


DELTA = os.environ.get("GEN_DELTA", "1") == "1"
MOVE0 = os.environ.get("GEN_MOVE0", "0") == "1"     # joint-sum planner joint column 0: S1(e)[0] skips T1, enters a' and e'
MOVED = []
J12 = set(int(v) for v in os.environ.get("GEN_J12", "").split(",") if v)   # rounds built by joint-sum planner's jmadd.build2
os.environ.setdefault("JOINT_LAMBDA", os.environ.get("GEN_LAMBDA", "0"))
import jmadd
PCJ = {}
FUSEF = os.environ.get("GEN_FUSEF", "1") == "1"          # F-chains inside the fused rounds 62/63
TAILPLAIN = os.environ.get("GEN_TAILPLAIN", "0") == "1"  # rounds 62/63 use the old Maj cells (easier tail lemma)


def maj_item(x, a, b, c):
    """rooted Maj: (roots a, raws M) with Maj = a ^ M, M = (a^b) & ~(b^c); None if anything is constant."""
    if not DELTA or any(l <= 1 for w in (a, b, c) for l in w): return None
    Xf = lambda p, q: gen3.X(x, p, q)
    M = [gen3.A(x, madd2.fdiff(Xf, a[i], b[i]), madd2.fdiff(Xf, b[i], c[i]) ^ 1) for i in range(32)]
    return (a, M) if all(m > 1 for m in M) else None


def choose(t1rows, t1pairs, t1f, Kt, s0, mj, d, rng, iters=60, ka=0, ke=0, ev=None):
    """ev: optional extra row (bit 0 only) moved out of T1 into both a' and e' (joint-sum planner joint column 0)."""
    V = [2] * 32
    X0 = [[2] + [0] * 31] if ev is not None else []
    NEG = list(itertools.product((False, True), repeat=3))
    def rc(C, n):
        ea = E([V, s0] + X0, C + ka, (), n[2], (), False, True) if mj is None else E([V, s0, mj] + X0, C + ka, (), n[2])
        return (E(t1rows, Kt - C, t1pairs, n[0], t1f) + ea + E([d, V] + X0, C + ke, (), n[1]))
    def best_neg(C): return min((rc(C, n), n) for n in NEG)
    bc = 0; bv, bn = best_neg(0); base = bv
    for it in range(iters):
        C = rng.getrandbits(32) if it % 2 else bc ^ (1 << rng.randrange(32))
        v, n = best_neg(C)
        if v < bv: bv, bc, bn = v, C, n
    imp = True
    while imp:
        imp = False
        for i in range(32):
            for C in (bc ^ (1 << i), (bc + (1 << i)) & M32, (bc - (1 << i)) & M32):
                v, n = best_neg(C)
                if v < bv: bv, bc, bn, imp = v, C, n, True
    return base, bv, bc, bn


JS = os.environ.get("GEN_JS", "")       # joint sigma0/sigma1 SLP file (lines "out a b", masks over the word bits)
_slp = [tuple(int(v) for v in l.split()) for l in open(JS)] if JS else None
_sig = {}


def joint_sigma(x, b, W, t):
    """(sigma0(W_t), sigma1(W_t)) for a word that needs both (linear-layer optimiser's joint program), else the plain builders."""
    if t in _sig: return _sig[t]
    if _slp is None or not (14 <= t <= 48):
        _sig[t] = (b.s0(W[t]) if 1 <= t <= 48 else None, b.s1(W[t]) if 14 <= t <= 61 else None)
        return _sig[t]
    lit = {1 << i: W[t][i] for i in range(32)}
    for o, p, q in _slp: lit[o] = gen3.X(x, lit[p], lit[q])
    def tap(k, rots, sh):
        m = 0
        for r in rots: m |= 1 << ((k + r) % 32)
        if k + sh < 32: m |= 1 << (k + sh)
        return m
    _sig[t] = ([lit[tap(k, (7, 18), 3)] for k in range(32)], [lit[tap(k, (17, 19), 10)] for k in range(32)])
    return _sig[t]


def cover_set():
    """27 words covering every schedule sum t = 16..63 (t-7 or t-16 in J)."""
    J = set()
    for t0 in range(16, 25):
        chain = list(range(t0, 64, 9))
        i = 0
        while i < len(chain):
            if i + 1 < len(chain): J.add(chain[i] - 7); i += 2      # word t-7 = (t+9)-16 covers t and t+9
            else: J.add(chain[i] - 7 if chain[i] - 7 <= 56 and chain[i] - 7 >= 9 else chain[i] - 16); i += 1
    return J


ANAT = []


def build(J, log=None):
    rng = random.Random(20261118)
    ANAT.clear()
    _last = [0, 0]
    def mark(tag):
        na = sum(1 for nd in x.nodes if nd[0] == "AND"); ng = len(x.nodes)
        ANAT.append((tag, na - _last[0], ng - _last[1])); _last[0], _last[1] = na, ng
    _sig.clear(); madd2.FREG.clear()
    x = XAG(512); b = FB(x)
    W = g2.inputs_words(x, 16)
    mark("in")
    Wrows = {}
    for t in range(16, 64):
        s1, s0 = joint_sigma(x, b, W, t - 2)[1], joint_sigma(x, b, W, t - 15)[0]
        Wrows[t] = [s1, W[t-7], s0, W[t-16]]
        if t >= 62:
            W.append(None); continue
        opts = [((), [s1, W[t-7], s0, W[t-16]])]
        if t - 7 in J: opts.append(((W[t-7],), [s1, s0, W[t-16]]))
        if t - 16 in J: opts.append(((W[t-16],), [s1, W[t-7], s0]))
        if t - 7 in J and t - 16 in J: opts.append(((W[t-7], W[t-16]), [s1, s0]))
        wf = t in J
        fr, rows = min(opts, key=lambda o: E(o[1], 0, (), False, o[0], wf))
        W.append(S(x, rows, 0, (), False, fr, wf))
        mark("W%d" % t)
    st = [g2.const_word(v) for v in sref.IV]
    for t in range(64):
        a, bb, c, d, e, f, gg, h = st
        fuse = W[t] is None
        if fuse:
            fw = [W[j] for j in (t - 7, t - 16) if j in J and T1F and FUSEF]
            wr = [Wrows[t][0], Wrows[t][2]] + [W[j] for j in (t - 7, t - 16) if not (j in J and T1F and FUSEF)]
            t1f = tuple(fw)
        else:
            t1f = (W[t],) if (t in J and T1F) else ()
            wr = [] if t1f else [W[t]]
        t1rows, t1pairs = [b.S1(e)] + wr, [(h, chpair(x, e, f, gg, h))]
        last = (t == 63)
        ka = sref.IV[0] if last else 0
        ke = sref.IV[4] if last else 0
        Kt = sref.K[t]
        s0 = b.S0(a)
        if t in J12:      # joint-sum planner joint (a', e') column compression, explicit Maj, no F-chains / self-delta
            an, en = jmadd.build2(x, t1rows + list(t1f), t1pairs, [s0, b.maj(a, bb, c)], [d], Kt, ka, ke, cache=PCJ)
            st = [an, a, bb, c, en, e, f, gg]
            mark("R%d" % t)
            continue
        dw = maj_item(x, a, bb, c) if (t < 62 or not TAILPLAIN) else None
        mj = None if dw else b.maj(a, bb, c)
        if dw is None:
            sep = E(t1rows, Kt, t1pairs, False, t1f) + E([[2] * 32, s0, mj], ka) + E([d, [2] * 32], ke)
        else:
            sep = 0; inl = 1
        allrows = t1rows + [w for w in t1f]
        if dw is None:
            inl = E(allrows + [s0, mj], Kt + ka, t1pairs) + E(allrows + [d], Kt + ke, t1pairs)
        if inl < sep:
            an = S(x, allrows + [s0, mj], Kt + ka, t1pairs)
            en = S(x, allrows + [d], Kt + ke, t1pairs)
        else:
            cands = [(choose(t1rows, t1pairs, t1f, Kt, s0, mj, d, rng, ka=ka, ke=ke), t1rows, [])]
            if MOVE0 and t1rows[0][0] > 1:
                r0 = [0] + list(t1rows[0][1:]); evr = [t1rows[0][0]] + [0] * 31
                cands.append((choose([r0] + t1rows[1:], t1pairs, t1f, Kt, s0, mj, d, rng, ka=ka, ke=ke, ev=evr),
                              [r0] + t1rows[1:], [evr]))
            (base, bv, bc, bn), t1r, X0 = min(cands, key=lambda z: z[0][1])
            if log is not None: log.append((t, base, bv, hex(bc), bn, bool(X0)))
            MOVED.append(bool(X0))
            t1 = S(x, t1r, Kt - bc, t1pairs, bn[0], t1f)
            an = S(x, [t1, s0, mj] + X0, bc + ka, (), bn[2]) if dw is None else S(x, [t1, s0] + X0, bc + ka, (), bn[2], (), False, dw)
            en = S(x, [d, t1] + X0, bc + ke, (), bn[1])
        st = [an, a, bb, c, en, e, f, gg]
        mark("R%d" % t)
    out = [st[0]] + [S(x, [st[k]], sref.IV[k]) for k in (1, 2, 3)] + [st[4]] + \
          [S(x, [st[k]], sref.IV[k]) for k in (5, 6, 7)]
    mark("final")
    x.outputs = [w[31 - j] for w in out for j in range(32)]
    return x


def parse_J(s):
    if s == "all": return set(range(0, 62))
    if s == "none": return set()
    if s == "cover": return cover_set()
    if s.startswith("cover+"): return cover_set() | set(int(v) for v in s[6:].split(",") if v)
    return set(int(v) for v in s.split(",") if v)


if __name__ == "__main__":
    J = parse_J(sys.argv[2] if len(sys.argv) > 2 else "cover")
    x = build(J)
    save_cache()
    y = cost.canonicalize(x)
    r = cost.report_xag(y); r["J"] = sorted(J)
    print(json.dumps(r))
    if os.environ.get("GEN_ANAT"): json.dump(ANAT, open(os.environ["GEN_ANAT"], "w"))
    if len(sys.argv) > 1 and sys.argv[1] != "-": open(sys.argv[1], "w").write(write_xag(y))
