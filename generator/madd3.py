"""madd3: carry-root multi-operand adder planner (madd2 design) with the joint-sum planner AND-aware objective (cost = gates + LAM * ANDs, env GEN_LAMBDA).

carry-root generator madd2: column-compression generator/constant-redistribution generator madd column compression + *rooted raw carries* ("F-chains").

New idea (carry-root generator).  A full-adder carry can be written  c = u ^ N  with N one AND gate and u ("root") any of the
three inputs.  If the next column pairs c with a bit q, the pair's second component is
    c ^ q = N ^ (u ^ q),
so when the form u ^ q is available the explicit carry XOR and the pairing XOR collapse into ONE XOR.
For a plain operand word X take u = X[i-1], q = X[i]:  F_i = X[i-1] ^ X[i] (the "difference word" of X).
F depends only on X, so a word that is an operand of several sums pays for F once and every sum that
runs an F-chain along X saves one gate per column (4 gates per elimination with no pairing cost).

Column items: singles S, pairs P (u, w), constant k, and per F-row f: its plain bit x_f and possibly an
incoming raw carry r_f (the AND output N of column i-1; the carry is X_f[i-1] ^ N).
Extra ops (cost in gates, F itself not counted):
  X2S   x_f -> s                               0
  RAW2S r_f -> s            (c = N ^ X[i-1])    1
  PAIRF x_f, r_f -> P       (w = N ^ F_i)       1
  FAPF  x_f, r_f, s -> s + raw                 4
  FA3F  x_f, s, s   -> s + raw                 4
  HAF   x_f, s      -> s + raw                 2
  HAPF  x_f, r_f    -> s + raw                 2
  FAKF  x_f, s, 1   -> s + raw                 2
  FAPFK x_f, r_f, 1 -> s + raw                 2

Second idea (carry-root generator): the difference word of a SUM OUTPUT is a free by-product ("self-F").  Root the carry of the
last elimination of column i-1 at the column output OUT[i-1] (maj = sum ^ (A | B), A, B two plane forms):
c = OUT[i-1] ^ N.  In column i let the last elimination be (s; u, c).  Then
    A = s ^ u,   F_i = A ^ N  (= OUT[i-1] ^ OUT[i]),   OUT[i] = F_i ^ OUT[i-1],   B = OUT[i] ^ s (= u ^ c),
    N_i = A | B                       -- 5 gates, no explicit carry, no pairing XOR, and F_i for free.
  SELFF  s, s, oraw -> OUT + oraw (+F_i)       5   (terminal op of a column)
  SELF0  s, s, s    -> OUT + oraw              4   (terminal; chain start)
  SELFH  s, s       -> OUT + oraw              2   (terminal)
  ORAW2S oraw -> s                             1

Third idea (carry-root generator): the time-difference word a' ^ a is a free by-product ("self-delta").  Maj(a,b,c) =
a ^ M with M = (a^b) & ~(b^c) one AND.  Give the a' column the wire Maj as the rooted item D = (a_i, M_i)
and let it enter the LAST elimination (s; u, Maj):
    A = s^u,  delta = A ^ M (= a'_i ^ a_i),  OUT = delta ^ a_i,  B = OUT ^ u (= s ^ Maj),  N = A & B
(carry = s ^ N).  The explicit Maj XOR disappears into OUT and delta is exactly the form a'^a that the next
round's Maj needs (so far 32 XOR per round).
  MDFAD s, s, P, D -> OUT + pair (+delta)     10   (terminal; = MDFA + PAIR + Maj-XOR, delta free)
  FAD   s, s, D    -> OUT + carry (+delta)      6   (terminal)
  HAD   s, D       -> OUT + carry (+delta)      3   (terminal)
  FAKD  s, D, 1    -> OUT + carry (+delta)      3   (terminal)
  D0    D [,1]     -> OUT [+ carry] (+delta)    1   (terminal)
  D2S   D -> s  (explicit Maj)                  1   (delta then costs 1 elsewhere: planner penalty)
"""
import sys, os, functools
HERE = os.path.dirname(os.path.abspath(__file__))
import madd
from madd import OPS, CELLS, EXTRA, Cells

N = 32
SMAX, PMAX = 5, 4
FOPS = {  # name: (dS, dP, needk, cost, needS, need_raw, emits_raw)
    "X2S": (1, 0, 0, 0, 0, 0, 0),
    "RAW2S": (1, 0, 0, 1, 0, 1, 0),      # does not consume x
    "PAIRF": (0, 1, 0, 1, 0, 1, 0),
    "FAPF": (0, 0, 0, 4, 1, 1, 1),
    "FA3F": (-1, 0, 0, 4, 2, 0, 1),
    "HAF": (0, 0, 0, 2, 1, 0, 1),
    "HAPF": (1, 0, 0, 2, 0, 1, 1),
    "FAKF": (0, 0, 1, 2, 1, 0, 1),
    "FAPFK": (1, 0, 1, 2, 0, 1, 1),
}
ENABLED = set(os.environ.get("MADD2_OPS", ",".join(FOPS)).split(","))
SELF = os.environ.get("MADD2_SELF", "1") == "1"
DELTA = os.environ.get("MADD2_DELTA", "1") == "1"
USEREG = os.environ.get("MADD2_REG", "1") == "1"   # 0: consumers rebuild X[i-1]^X[i] explicitly (proof variant B)
FREG = {}
LAM = float(os.environ.get("GEN_LAMBDA", "0"))
ANDS = {"MDFA": 2, "MDFAK": 2, "FA3": 1, "FAP": 1, "FAK": 1, "FAPK": 1, "HA": 1, "HAP": 1, "HAK": 0,
        "PAIR": 0, "MDFA0": 3, "MDFAZ": 3, "FA3KP": 1, "FAPKP": 1, "Q4": 3, "MDFAKS": 3,
        "X2S": 0, "RAW2S": 0, "PAIRF": 0, "FAPF": 1, "FA3F": 1, "HAF": 1, "HAPF": 1, "FAKF": 1, "FAPFK": 1,
        "SELFF": 1, "SELF0": 1, "SELFH": 1, "ORAW2S": 0, "MDFAD": 2, "FAD": 1, "HAD": 1, "FAKD": 1, "D0": 0,
        "D2S": 0}


def W(name, c): return c + LAM * ANDS[name]      # (lit_a & ~1, lit_b & ~1) sorted -> literal of (a & ~1) ^ (b & ~1) when it exists as a by-product


def fdiff(X, a, b):
    """literal of a ^ b, using a registered self-F by-product when there is one."""
    key = tuple(sorted((a & ~1, b & ~1)))
    if USEREG and key in FREG: return FREG[key] ^ ((a ^ b) & 1)
    return X(a, b)


@functools.lru_cache(maxsize=None)
def col_options(S, P, k, xm, rm, o, d, last):
    """dict (emit_s, emit_p, rawout, oraw_out, fdone, ddone) -> (cost, ops).
    xm/rm: bitmasks of unconsumed plain bits / raw carries; o: unconsumed output-rooted raw carry;
    d: unconsumed rooted item D (Maj)."""
    res = {}
    nx = bin(xm).count("1")
    def put(key, tot, ops):
        if key not in res or tot < res[key][0]: res[key] = (tot, ops)
    if P == 0 and rm == 0 and not o and not d and S + k + nx <= 1:
        put((0, 0, 0, 0, 0, 0), 0, ())
    if SELF and not d and P == 0 and rm == 0 and k == 0 and not last:
        if o and S + nx == 2: put((0, 0, 0, 1, 1, 0), W("SELFF", 5), ("SELFF",))
        if not o and S + nx == 3: put((0, 0, 0, 1, 0, 0), W("SELF0", 4), ("SELF0",))
        if not o and S + nx == 2: put((0, 0, 0, 1, 0, 0), W("SELFH", 2), ("SELFH",))
    if d and DELTA and rm == 0 and not o and not last:
        if S + nx == 2 and P == 1 and k == 0: put((0, 1, 0, 0, 0, 1), W("MDFAD", 10), ("MDFAD",))
        if S + nx == 2 and P == 0 and k == 0: put((1, 0, 0, 0, 0, 1), W("FAD", 6), ("FAD",))
        if S + nx == 1 and P == 0 and k == 0: put((1, 0, 0, 0, 0, 1), W("HAD", 3), ("HAD",))
        if S + nx == 1 and P == 0 and k == 1: put((1, 0, 0, 0, 0, 1), W("FAKD", 3), ("FAKD",))
        if S + nx == 0 and P == 0: put((k, 0, 0, 0, 0, 1), 1, ("D0",))
    def upd(sub, c, es, ep, ro, name):
        for (a, b, r, oo, fd, dd), (cc, ops) in sub.items():
            if r & ro: continue
            put((a + es, b + ep, r | ro, oo, fd, dd), cc + c, (name,) + ops)
    for name, (dS, dP, nk, c, es, ep, nS, nP) in OPS.items():
        if S < nS or P < nP or (nk and not k): continue
        upd(col_options(S + dS, P + dP, k - nk, xm, rm, o, d, last), W(name, c), es, ep, 0, name)
    if o:
        upd(col_options(S + 1, P, k, xm, rm, 0, d, last), 1, 0, 0, 0, "ORAW2S")
    if d:
        upd(col_options(S + 1, P, k, xm, rm, o, 0, last), 1, 0, 0, 0, "D2S")
    f = 0
    while (xm | rm) >> f:
        bit = 1 << f
        for name, (dS, dP, nk, c, nS, nr, er) in FOPS.items():
            if name not in ENABLED: continue
            if nr and not rm & bit: continue
            if name != "RAW2S" and not xm & bit: continue
            if S < nS or (nk and not k): continue
            if er and last: continue
            nxm = xm if name == "RAW2S" else xm & ~bit
            nrm = rm & ~bit if nr else rm
            if name != "RAW2S" and not nr and rm & bit:
                continue          # consuming x_f while its raw carry is pending: force raw first (RAW2S/PAIRF)
            upd(col_options(S + dS, P + dP, k - nk, nxm, nrm, o, d, last), W(name, c), 0, 0, bit if er else 0, (name, f))
        f += 1
    return res


def top_cost(S, P, k, xm, rm, o, d):
    """last column: XOR of everything.  (#items - 1) + (#raw carries) + 2 if an output-rooted raw carry comes
    in; a rooted item D counts as two items (M and its root; the root is XORed last, delta is free)."""
    nx = bin(xm).count("1"); nr = bin(rm).count("1")
    items = S + P + nx + (2 if d else 0)
    return max(0, items - 1) + nr + (2 if o else 0)


def plan(raw, kbits, rawp, nf, n=N, wantF=False, dcols=None):
    full = (1 << nf) - 1
    dcols = dcols or [0] * n
    states = {(0, 0, 0, 0): (0, [])}
    for i in range(n):
        new = {}
        for (s, p, r, o), (c, hist) in states.items():
            S, P, k = raw[i] + s, p + rawp[i], kbits[i]
            if i == n - 1:
                tot = c + top_cost(S, P, k, full, r, o, dcols[i]) + (1 if wantF and not o else 0)
                if (0, 0, 0, 0) not in new or tot < new[(0, 0, 0, 0)][0]:
                    new[(0, 0, 0, 0)] = (tot, hist + [("TOP",)])
                continue
            for (es, ep, ro, oo, fd, dd), (cc, ops) in col_options(S, P, k, full, r, o, dcols[i], False).items():
                if es > SMAX or ep > PMAX: continue
                tot = c + cc + (1 if wantF and i >= 1 and not fd else 0) + (1 if dcols[i] and not dd else 0)
                key = (es, ep, ro, oo)
                if key not in new or tot < new[key][0]:
                    new[key] = (tot, hist + [ops])
        states = new
    return states[(0, 0, 0, 0)]


def _key(rows, K, n, pairs, frows, wantF, dcols=None):
    cols, pcols, K = madd.prep(rows, K, n, pairs)
    raw = [len(cols[i]) for i in range(n)]
    rawp = [len(pcols[i]) for i in range(n)]
    kb = [(K >> i) & 1 for i in range(n)]
    return cols, pcols, kb, (tuple(raw), tuple(kb), tuple(rawp), len(frows), bool(wantF), SELF, tuple(dcols) if dcols else None, LAM)


def estimate(rows, K=0, n=N, pairs=(), frows=(), plan_cache=None, wantF=False, dcols=None):
    if not frows and not wantF and not dcols and LAM == 0:
        return madd.estimate(rows, K, n, pairs, plan_cache)
    cols, pcols, kb, key = _key(rows, K, n, pairs, frows, wantF, dcols)
    if plan_cache is not None and key in plan_cache: return plan_cache[key][0]
    r = plan(list(key[0]), kb, list(key[2]), len(frows), n, wantF, dcols)
    if plan_cache is not None: plan_cache[key] = r
    return r[0]


def build_sum(x, rows, K=0, n=N, plan_cache=None, pairs=(), frows=(), XF=None, AF=None, wantF=False, dwire=None):
    """frows: list of plain words (n variable literals each).  XF/AF: folding XOR/AND builders (x, a, b).
    wantF: the difference word of the output will be used; it is registered in FREG."""
    if not frows and not wantF and not dwire and LAM == 0:
        return madd.build_sum(x, rows, K, n, plan_cache, pairs)
    for fw in frows: assert all(l > 1 for l in fw)
    dcols = [1] * n if dwire else None      # dwire = (roots, raws): Maj[i] = roots[i] ^ raws[i], all variable
    c = Cells(x)
    X = (lambda a, b: XF(x, a, b)) if XF else c.X
    A = (lambda a, b: AF(x, a, b)) if AF else c.A
    cols, pcols, kb, key = _key(rows, K, n, pairs, frows, wantF, dcols)
    if plan_cache is not None and key in plan_cache:
        cost, ops = plan_cache[key]
    else:
        cost, ops = plan(list(key[0]), kb, list(key[2]), len(frows), n, wantF, dcols)
        if plan_cache is not None: plan_cache[key] = (cost, ops)
    out = []
    cs, cp = [], []
    rawin = [None] * len(frows)
    oin = None
    def reg(prev, cur, F):
        FREG[tuple(sorted((prev & ~1, cur & ~1)))] = F ^ ((prev ^ cur) & 1)
    for i in range(n):
        S = list(cols[i]) + cs
        P = list(pcols[i]) + list(cp)
        k = kb[i]
        es, ep = [], []
        rawout = [None] * len(frows)
        oout = None
        xb = [fw[i] for fw in frows]
        def Fw(f):   # w = N ^ F_i  = carry ^ X[i]
            return X(rawin[f], fdiff(X, frows[f][i - 1], frows[f][i]))
        if ops[i] == ("TOP",):
            items = list(S) + [w for (u, w) in P]
            for f in range(len(frows)):
                items.append(Fw(f) if rawin[f] is not None else xb[f])
            acc = k
            for l in items: acc = X(acc, l)
            if oin is not None:
                F = X(acc, oin); acc = X(F, out[i - 1]); reg(out[i - 1], acc, F)
            if dwire:
                dl = X(acc, dwire[1][i]); acc = X(dl, dwire[0][i]); reg(dwire[0][i], acc, dl)
            out.append(acc)
            break
        done = False
        dpend = bool(dwire)
        for op in ops[i]:
            if op in ("MDFAD", "FAD", "HAD", "FAKD", "D0"):
                rest = S + [u for u in xb if u is not None]
                ra, M = dwire[0][i], dwire[1][i]
                assert all(r is None for r in rawin) and oin is None
                if op == "MDFAD":
                    x_, u = rest; (up, wp) = P.pop(0); assert not P and k == 0
                    s1 = X(x_, wp); g = X(x_, up); N1 = A(g ^ 1, wp ^ 1) ^ 1
                    Af = X(s1, u); dl = X(Af, M); o_ = X(dl, ra); B = X(o_, u); N2 = A(Af, B)
                    ep.append((X(s1, N1), X(N1, N2)))
                elif op == "FAD":
                    s_, u = rest; assert not P and k == 0
                    Af = X(s_, u); dl = X(Af, M); o_ = X(dl, ra); B = X(o_, u); es.append(X(s_, A(Af, B)))
                elif op == "HAD":
                    (s_,) = rest; assert not P and k == 0
                    dl = X(s_, M); o_ = X(dl, ra); es.append(A(s_, o_ ^ 1))
                elif op == "FAKD":
                    (s_,) = rest; assert not P and k == 1
                    dl0 = X(s_, M); on = X(dl0, ra); o_ = on ^ 1; dl = dl0 ^ 1; es.append(A(s_ ^ 1, on ^ 1) ^ 1); k = 0
                else:
                    assert not rest and not P
                    mj_ = X(M, ra); o_ = mj_ ^ k; dl = M ^ k
                    if k: es.append(mj_)
                    k = 0
                reg(ra, o_, dl)
                out.append(o_); done = True; dpend = False
                break
            if op == "D2S":
                S.append(X(dwire[0][i], dwire[1][i])); dpend = False
                continue
            if op in ("SELFF", "SELF0", "SELFH"):
                rest = S + [u for u in xb if u is not None]
                assert not P and all(r is None for r in rawin) and k == 0
                if op == "SELFF":
                    s_, u = rest; Af = X(s_, u); F = X(Af, oin); o_ = X(F, out[i - 1]); reg(out[i - 1], o_, F)
                    B = X(o_, s_); oout = A(Af ^ 1, B ^ 1) ^ 1; oin = None
                elif op == "SELF0":
                    s_, u, v = rest; Af = X(s_, u); o_ = X(Af, v); B = X(o_, s_); oout = A(Af ^ 1, B ^ 1) ^ 1
                else:
                    s_, u = rest; o_ = X(s_, u); oout = A(s_ ^ 1, u ^ 1) ^ 1
                out.append(o_); done = True
                break
            if op == "ORAW2S":
                S.append(X(oin, out[i - 1])); oin = None
            elif isinstance(op, tuple):
                name, f = op
                u = xb[f]
                if name == "X2S": S.append(u)
                elif name == "RAW2S": S.append(X(rawin[f], frows[f][i - 1])); rawin[f] = None
                elif name == "PAIRF": P.append((u, Fw(f))); rawin[f] = None
                elif name == "FAPF":
                    w = Fw(f); rawin[f] = None; a = S.pop(0)
                    S.append(X(w, a)); rawout[f] = A(w, X(u, a))
                elif name == "FA3F":
                    a, b = S.pop(0), S.pop(0); t1 = X(u, a); t2 = X(u, b)
                    S.append(X(t1, b)); rawout[f] = A(t1, t2)
                elif name == "HAF":
                    a = S.pop(0); S.append(X(u, a)); rawout[f] = A(u, a ^ 1)
                elif name == "HAPF":
                    w = Fw(f); rawin[f] = None; S.append(w); rawout[f] = A(u, w)
                elif name == "FAKF":
                    a = S.pop(0); k = 0; S.append(X(u, a) ^ 1); rawout[f] = A(u ^ 1, a)
                elif name == "FAPFK":
                    w = Fw(f); rawin[f] = None; k = 0; S.append(w ^ 1); rawout[f] = A(u ^ 1, w)
                else: raise ValueError(op)
                if name != "RAW2S": xb[f] = None
            elif op == "MDFA":
                a = S.pop(0); p1 = P.pop(0); p2 = P.pop(0)
                sm, pr = c.mdfa(a, p1, p2); S.append(sm); ep.append(pr)
            elif op == "MDFAK":
                p1 = P.pop(0); p2 = P.pop(0)
                sm, pr = c.mdfa(1, p1, p2); S.append(sm); ep.append(pr); k = 0
            elif op == "FA3":
                a, b, d = S.pop(0), S.pop(0), S.pop(0)
                t = c.X(a, d); sm = c.X(t, b); cy = c.X(c.A(t, c.X(b, d)), d)
                S.append(sm); es.append(cy)
            elif op == "FAP":
                (u, w) = P.pop(0); a = S.pop(0)
                sm = c.X(w, a); cy = c.X(u, c.A(w, c.X(u, a)))
                S.append(sm); es.append(cy)
            elif op == "FAK":
                a, b = S.pop(0), S.pop(0); k = 0
                S.append(c.X(a, b) ^ 1); es.append(c.A(a ^ 1, b ^ 1) ^ 1)
            elif op == "FAPK":
                (u, w) = P.pop(0); k = 0
                S.append(w ^ 1); es.append(c.A(u ^ 1, w ^ 1) ^ 1)
            elif op == "HA":
                a, b = S.pop(0), S.pop(0)
                S.append(c.X(a, b)); es.append(c.A(a, b))
            elif op == "HAP":
                (u, w) = P.pop(0)
                S.append(w); es.append(c.A(u, w ^ 1))
            elif op == "HAK":
                a = S.pop(0); k = 0
                S.append(a ^ 1); es.append(a)
            elif op in EXTRA:
                cl = CELLS[op]; wires = []
                for t in cl["spec"].split(","):
                    if t == "s": wires.append(S.pop(0))
                    elif t == "p": wires += list(P.pop(0))
                    else: k = 0
                o = c.cell(cl, wires); S.append(o[0]); ep.append((o[1], o[2]))
                if len(o) > 3: es.append(o[3])
            elif op == "PAIR":
                a, b = S.pop(0), S.pop(0)
                P.append((a, c.X(a, b)))
            else:
                raise ValueError(op)
        assert len(P) == 0 and all(r is None for r in rawin) and oin is None and not dpend, (i, S, P, k, rawin)
        if not done:
            rest = [u for u in xb if u is not None]
            assert len(S) + len(rest) + k <= 1, (i, S, rest, k)
            out.append(S[0] if S else (rest[0] if rest else k))
        cs, cp, rawin, oin = es, ep, rawout, oout
    return out


if __name__ == "__main__":
    for nm, raw, rawp, nf, wf in (("sched nf=0", 4, 0, 0, 0), ("sched nf=1", 3, 0, 1, 0), ("sched nf=2", 2, 0, 2, 0),
                                  ("sched nf=0 wantF", 4, 0, 0, 1), ("sched nf=1 wantF", 3, 0, 1, 1), ("sched nf=2 wantF", 2, 0, 2, 1),
                                  ("2op nf=0", 2, 0, 0, 0), ("2op nf=0 wantF", 2, 0, 0, 1), ("2op nf=1", 1, 0, 1, 0), ("2op nf=1 wantF", 1, 0, 1, 1),
                                  ("3op nf=0", 3, 0, 0, 0), ("3op nf=0 wantF", 3, 0, 0, 1), ("3op nf=1", 2, 0, 1, 0), ("3op nf=1 wantF", 2, 0, 1, 1),
                                  ("T1 nf=0", 2, 1, 0, 0), ("T1 nf=1", 1, 1, 1, 0), ("T1 nf=1 wantF", 1, 1, 1, 1)):
        print(nm, plan([raw] * N, [0] * N, [rawp] * N, nf, N, bool(wf))[0])
    print("a' with explicit Maj (3 singles) + 31 Maj XOR + 32 delta:", plan([3] * N, [0] * N, [0] * N, 0)[0] + 63)
    print("a' with D item (2 singles + D):", plan([2] * N, [0] * N, [0] * N, 0, N, False, [1] * N)[0])
