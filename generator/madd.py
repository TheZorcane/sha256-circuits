"""column-compression generator: multi-operand addition mod 2^32 by column compression with MDFA cells and
encoded carry pairs, planned by an exact DP over a fixed op set.

Column items: singles (literals), pairs (u, w) standing for the two bits u and u^w
(same weight), and at most one constant-1 bit (constants of all rows folded into K).

Ops inside a column (cost = gates before strash/folding):
  MDFA  (s; P; P) -> s + pair to next col           8   [GKL arXiv:2509.13966 / cell study exact]
  MDFAK (1; P; P) -> s + pair                       6   (MDFA with x3 = 1, constant folded)
  FA3   (s,s,s)   -> s + carry                      5
  FAP   (P; s)    -> s + carry                      4   sum = w^s, carry = u ^ (w & (u^s))
  FAK   (s,s,1)   -> s + carry                      2   sum = ~(x^y), carry = x|y
  FAPK  (P; 1)    -> s + carry                      1   sum = ~w, carry = u|w
  HA    (s,s)     -> s + carry                      2
  HAP   (P)       -> s + carry                      1   sum = w, carry = u & ~w
  HAK   (s; 1)    -> s + carry                      0   sum = ~x, carry = x
  PAIR  (s,s)     -> P                              1
Last column (31): XOR of all singles and pair second components.
"""
import sys, os, functools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))

N = 32
OPS = {  # name: (dS, dP, needk, cost, emit_s, emit_p, needS, needP)
    "MDFA": (0, -2, 0, 8, 0, 1, 1, 2),
    "MDFAK": (1, -2, 1, 6, 0, 1, 0, 2),
    "FA3": (-2, 0, 0, 5, 1, 0, 3, 0),
    "FAP": (0, -1, 0, 4, 1, 0, 1, 1),
    "FAK": (-1, 0, 1, 2, 1, 0, 2, 0),
    "FAPK": (1, -1, 1, 1, 1, 0, 0, 1),
    "HA": (-1, 0, 0, 2, 1, 0, 2, 0),
    "HAP": (1, -1, 0, 1, 1, 0, 0, 1),
    "HAK": (0, 0, 1, 0, 1, 0, 1, 0),
    "PAIR": (-2, 1, 0, 1, 0, 0, 2, 0),
}
# exact-synthesized cells (cells.json, column-compression generator cells_exact.py): spec -> s + carry pair
EXTRA = {
    "MDFA0": (-1, -1, 0, 7, 0, 1, 2, 1),
    "MDFAZ": (1, -2, 0, 6, 0, 1, 0, 2),
    "FA3KP": (-2, 0, 1, 4, 0, 1, 3, 0),
    "FAPKP": (0, -1, 1, 3, 0, 1, 1, 1),
    "Q4": (-3, 0, 0, 8, 0, 1, 4, 0),
    "MDFAKS": (0, -2, 1, 7, 1, 1, 1, 2),   # s,p,p,1 -> s + pair + single (E8)
}
import json as _json
CELLS = _json.load(open(os.path.join(HERE, "cells.json")))
USE_EXTRA = os.environ.get("MADD_EXTRA", "1") == "1"
if USE_EXTRA:
    OPS.update(EXTRA)


@functools.lru_cache(maxsize=None)
def col_options(S, P, k):
    """dict (emit_s, emit_p) -> (cost, ops) for reducing a column to <= 1 bit."""
    res = {}
    bits = S + 2 * P + k
    if bits <= 1:
        res[(0, 0)] = (0, ())
        return res
    for name, (dS, dP, nk, c, es, ep, nS, nP) in OPS.items():
        if S < nS or P < nP or (nk and not k):
            continue
        sub = col_options(S + dS, P + dP, k - nk)
        for (a, b), (cc, ops) in sub.items():
            key = (a + es, b + ep)
            tot = cc + c
            if key not in res or tot < res[key][0]:
                res[key] = (tot, (name,) + ops)
    return res


SMAX, PMAX = 5, 4
FREEPAIR = os.environ.get("MADD_FREEPAIR") == "1"
STATS = {"free": 0, "pair": 0}


def top_cost(S, P, k):
    return max(0, S + P - 1)


def plan(raw, kbits, n=N, rawp=None):
    """raw[i] = #variable single bits in column i, kbits[i] in {0,1}.
    Returns (cost, list of op tuples per column)."""
    rawp = rawp or [0] * n
    states = {(0, 0): (0, [])}
    for i in range(n):
        new = {}
        for (s, p), (c, hist) in states.items():
            S, P, k = raw[i] + s, p + rawp[i], kbits[i]
            if i == n - 1:
                tot = c + top_cost(S, P, k)
                key = (0, 0)
                if key not in new or tot < new[key][0]:
                    new[key] = (tot, hist + [("TOP",)])
                continue
            for (es, ep), (cc, ops) in col_options(S, P, k).items():
                if es > SMAX or ep > PMAX:
                    continue
                tot = c + cc
                key = (es, ep)
                if key not in new or tot < new[key][0]:
                    new[key] = (tot, hist + [ops])
        states = new
    return states[(0, 0)]


class Cells:
    def __init__(s, x): s.x = x
    def X(s, a, b): return s.x.XOR(a, b)
    def A(s, a, b): return s.x.AND(a, b)

    def xors(s):
        d = s.x.__dict__.setdefault("_xidx", [0, set()])
        base = s.x.n_inputs + 1
        while d[0] < len(s.x.nodes):
            op, a, b = s.x.nodes[d[0]]
            if op == "XOR": d[1].add(tuple(sorted((a & ~1, b & ~1))))
            d[0] += 1
        return d[1]

    def cell(s, cl, wires):
        v = list(wires)
        for j, kk, f in cl["gates"]:
            p, q = v[j], v[kk]
            f = tuple(f)
            if f == (1, 1, 0): r = s.X(p, q)
            elif f == (0, 0, 1): r = s.A(p, q)
            elif f == (1, 1, 1): r = s.A(p ^ 1, q ^ 1) ^ 1
            elif f == (0, 1, 0): r = s.A(p, q ^ 1)
            elif f == (1, 0, 0): r = s.A(p ^ 1, q)
            elif f == (1, 0, 1): r = q
            elif f == (0, 1, 1): r = p
            else: raise ValueError(f)
            v.append(r)
        return [v[sg] ^ ng for sg, ng in cl["outs"]]

    def mdfa(s, a, p1, p2):
        b, p = p1; o1, o2 = p2
        g0 = s.X(a, b); g1 = s.X(a, p); g2 = s.X(o1, g1)
        g3 = s.A(o2 ^ 1, g2); g4 = s.X(g1, g3)
        g5 = s.A(p ^ 1, g0 ^ 1) ^ 1
        g6 = s.X(g3, g5); g7 = s.X(o2, g1)
        return g7, (g4, g6)


def prep(rows, K=0, n=N, pairs=()):
    M = (1 << n) - 1
    K &= M
    cols = [[] for _ in range(n)]
    for r in rows:
        for i, l in enumerate(r):
            if l in (0, 1):
                K = (K + (l << i)) & M
            else:
                cols[i].append(l)
    pcols = [[] for _ in range(n)]
    for (U, Wd) in pairs:          # pair word: bits U[i] and U[i]^Wd[i]
        for i in range(n):
            u, w = U[i], Wd[i]
            if u in (0, 1):
                K = (K + (u << i)) & M
                v = w ^ u
                if v in (0, 1): K = (K + (v << i)) & M
                else: cols[i].append(v)
            elif w in (0, 1):
                if w == 1: K = (K + (1 << i)) & M
                elif i + 1 < n: cols[i + 1].append(u)
            else:
                pcols[i].append((u, w))
    return cols, pcols, K


def estimate(rows, K=0, n=N, pairs=(), plan_cache=None):
    cols, pcols, K = prep(rows, K, n, pairs)
    raw = [len(cols[i]) for i in range(n)]
    rawp = [len(pcols[i]) for i in range(n)]
    kb = [(K >> i) & 1 for i in range(n)]
    key = (tuple(raw), tuple(kb), tuple(rawp))
    if plan_cache is not None and key in plan_cache:
        return plan_cache[key][0]
    r = plan(raw, kb, n, rawp)
    if plan_cache is not None: plan_cache[key] = r
    return r[0]


def build_sum(x, rows, K=0, n=N, plan_cache=None, pairs=()):
    """rows: list of words (lists of n literals, LSB first; lits 0/1 are constants);
    pairs: list of (U, Wd) words standing for rows U and U^Wd.  Returns the sum mod 2^n."""
    c = Cells(x)
    cols, pcols, K = prep(rows, K, n, pairs)
    raw = [len(cols[i]) for i in range(n)]
    rawp = [len(pcols[i]) for i in range(n)]
    kb = [(K >> i) & 1 for i in range(n)]
    key = (tuple(raw), tuple(kb), tuple(rawp))
    if plan_cache is not None and key in plan_cache:
        cost, ops = plan_cache[key]
    else:
        cost, ops = plan(raw, kb, n, rawp)
        if plan_cache is not None:
            plan_cache[key] = (cost, ops)
    out = []
    cs, cp = [], []   # carries into current column
    for i in range(n):
        S = list(cols[i]) + cs
        P = list(pcols[i]) + list(cp)
        k = kb[i]
        es, ep = [], []
        if ops[i] == ("TOP",):
            acc = k
            for l in S:
                acc = l if acc in (0, 1) and acc == 0 else (l ^ 1 if acc == 1 else c.X(acc, l))
            for (u, w) in P:
                acc = w if acc == 0 else (w ^ 1 if acc == 1 else c.X(acc, w))
            out.append(acc)
            break
        for op in ops[i]:
            if op == "MDFA":
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
                pick = None
                if FREEPAIR:
                    xs = c.xors()
                    for ia in range(len(S)):
                        for ib in range(ia + 1, len(S)):
                            kk = tuple(sorted((S[ia] & ~1, S[ib] & ~1)))
                            if kk in xs: pick = (ia, ib); break
                        if pick: break
                if pick:
                    STATS["free"] += 1
                    b = S.pop(pick[1]); a = S.pop(pick[0])
                else:
                    a, b = S.pop(0), S.pop(0)
                STATS["pair"] += 1
                P.append((a, c.X(a, b)))
            else:
                raise ValueError(op)
        assert len(P) == 0 and len(S) + k <= 1, (i, S, P, k)
        out.append(S[0] if S else k)
        cs, cp = es, ep
    return out


if __name__ == "__main__":
    for k in (2, 3, 4, 5):
        print(k, plan([k] * N, [0] * N)[0])
