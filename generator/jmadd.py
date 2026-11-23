"""joint-sum planner: JOINT column compression for two sums that share addends ("membership heap").

  A = (shared rows) + (A rows) + const        E = (shared rows) + (E rows) + const       (mod 2^32)

Every bit of a column is shared (addend of both sums), A-only or E-only.  A cell applied to shared bits
produces shared outputs (sum bit in the column, carries in the next one).  Shared bits that are left over
(the "remainder") are copied, for free, into the A column and the E column, where each sum is finished
privately.  The classical round (compute T1, then add) is the special case "remainder = 1 bit in every
column"; column-compression generator's "inline" is "remainder = everything".  The planner is an exact DP over columns on the state
(shared carries, A carries, E carries, constant borrows) and ALSO chooses the split of the constant between
the shared part (cS) and the private parts (K - cS) exactly (constant-redistribution generator's cshift, here without random search).

Cells and their gate costs are column-compression generator's (madd.OPS / cells.json).  Objective: gates + LAMBDA * ANDs
(env JOINT_LAMBDA, default 0), so LAMBDA > 0 trades XORs for ANDs.
"""
import sys, os, functools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import madd
from madd import OPS, CELLS, Cells, prep

N = 32
M32 = (1 << 32) - 1
LAM = float(os.environ.get("JOINT_LAMBDA", "0"))
JOINT = os.environ.get("JOINT_MODE", "1") == "1"        # 0: force full resolution of the shared part (T1 word)
BEAM = int(os.environ.get("JOINT_BEAM", "300"))
RMAX = int(os.environ.get("JOINT_RMAX", "6"))            # max remainder bits when not fully resolved
ANDS = {"MDFA": 2, "MDFAK": 2, "FA3": 1, "FAP": 1, "FAK": 1, "FAPK": 1, "HA": 1, "HAP": 1, "HAK": 0,
        "PAIR": 0, "MDFA0": 3, "MDFAZ": 3, "FA3KP": 1, "FAPKP": 1, "Q4": 3, "MDFAKS": 3}
SMAX, PMAX = 5, 4
POT = float(os.environ.get("JOINT_POT", "4")) + LAM


def wcost(op):
    return OPS[op][3] + LAM * ANDS[op]


@functools.lru_cache(maxsize=None)
def full_options(S, P, k):
    """(emit_s, emit_p) -> (cost, ops): reduce a column to <= 1 bit (madd.col_options, weighted cost)."""
    res = {}
    if S + 2 * P + k <= 1:
        res[(0, 0)] = (0, ())
        return res
    for name, (dS, dP, nk, c, es, ep, nS, nP) in OPS.items():
        if S < nS or P < nP or (nk and not k):
            continue
        for (a, b), (cc, ops) in full_options(S + dS, P + dP, k - nk).items():
            key = (a + es, b + ep)
            if key[0] > SMAX or key[1] > PMAX: continue
            tot = cc + wcost(name)
            if key not in res or tot < res[key][0]:
                res[key] = (tot, (name,) + ops)
    return res


@functools.lru_cache(maxsize=None)
def part_options(S, P, k):
    """(S', P', k', emit_s, emit_p) -> (cost, ops): apply any sequence of cells, stop anywhere."""
    res = {(S, P, k, 0, 0): (0, ())}
    for name, (dS, dP, nk, c, es, ep, nS, nP) in OPS.items():
        if S < nS or P < nP or (nk and not k):
            continue
        for (S1, P1, k1, a, b), (cc, ops) in part_options(S + dS, P + dP, k - nk).items():
            key = (S1, P1, k1, a + es, b + ep)
            if key[3] > SMAX or key[4] > PMAX: continue
            tot = cc + wcost(name)
            if key not in res or tot < res[key][0]:
                res[key] = (tot, (name,) + ops)
    return res


def f1(n): return max(0, n - 1)


@functools.lru_cache(maxsize=None)
def trans(cs, raw, kS, kA, kE, joint):
    """one column: carry state cs = (sS,pS,sA,pA,sE,pE), raw = variable bits (rS,rpS,rA,rpA,rE,rpE).
    returns {new carry state: (cost, (opsS, kpass, opsA, opsE))}"""
    sS, pS, sA, pA, sE, pE = cs
    rS, rpS, rA, rpA, rE, rpE = raw
    out = {}
    for (S1, P1, k1, es, ep), (c0, ops0) in part_options(rS + sS, rpS + pS, kS).items():
        if k1 and (kA or kE): continue
        rem = S1 + 2 * P1 + k1
        if rem > (RMAX if joint else 1): continue
        fa = full_options(S1 + rA + sA, P1 + rpA + pA, kA or k1)
        fe = full_options(S1 + rE + sE, P1 + rpE + pE, kE or k1)
        for (a_s, a_p), (ca, opsa) in fa.items():
            for (e_s, e_p), (ce, opse) in fe.items():
                key = (es, ep, a_s, a_p, e_s, e_p)
                tot = c0 + ca + ce
                if key not in out or tot < out[key][0]:
                    out[key] = (tot, (ops0, k1, opsa, opse))
    return out


def plan2(raws, KA, KE, n=N, joint=None, fixcS=None):
    """raws[i] = (rS,rpS,rA,rpA,rE,rpE) variable bit counts of column i; KA/KE total constants of the sums.
    Chooses cS (shared constant) bit by bit; private constants are KA - cS, KE - cS.
    returns (cost, [per column (kS, kA, kE, opsS, kpass, opsA, opsE)])"""
    if joint is None: joint = JOINT
    states = {((0,) * 6, 0, 0): (0, None)}
    back = []
    for i in range(n):
        new = {}
        last = (i == n - 1)
        for (cs, bA, bE), (cost, _) in states.items():
            for kS in (0, 1):
                if fixcS is not None and kS != (fixcS >> i) & 1: continue
                tA = ((KA >> i) & 1) - kS - bA; kA = tA & 1; nbA = 1 if tA < 0 else 0
                tE = ((KE >> i) & 1) - kS - bE; kE = tE & 1; nbE = 1 if tE < 0 else 0
                if last:
                    sS, pS, sA, pA, sE, pE = cs
                    rS, rpS, rA, rpA, rE, rpE = raws[i]
                    nS = rS + sS + rpS + pS
                    tot = cost + f1(nS) + f1(rA + sA + rpA + pA + min(nS, 1)) + f1(rE + sE + rpE + pE + min(nS, 1))
                    key = ((0,) * 6, 0, 0)
                    if key not in new or tot < new[key][0]:
                        new[key] = (tot, ((cs, bA, bE), (kS, kA, kE, None)))
                    continue
                for ncs, (cc, info) in trans(cs, raws[i], kS, kA, kE, joint).items():
                    key = (ncs, nbA, nbE)
                    tot = cost + cc
                    if key not in new or tot < new[key][0]:
                        new[key] = (tot, ((cs, bA, bE), (kS, kA, kE) + info))
        if len(new) > BEAM:
            # potential: every pending carry bit still has to be eliminated (about 4 gates + 1 AND each)
            keep = sorted(new.items(), key=lambda kv: kv[1][0] + POT * (
                kv[0][0][0] + kv[0][0][2] + kv[0][0][4] + 2 * (kv[0][0][1] + kv[0][0][3] + kv[0][0][5])))[:BEAM]
            new = dict(keep)
        back.append(new)
        states = new
    key = ((0,) * 6, 0, 0)
    cost = states[key][0]
    hist = []
    for i in range(n - 1, -1, -1):
        _, (prev, info) = back[i][key]
        hist.append(info)
        key = prev
    hist.reverse()
    return cost, hist


def plan1(raw, kb, n=N, rawp=None):
    """single sum, weighted objective (drop-in for madd.plan)."""
    rawp = rawp or [0] * n
    states = {(0, 0): (0, [])}
    for i in range(n):
        new = {}
        for (s, p), (c, hist) in states.items():
            S, P, k = raw[i] + s, p + rawp[i], kb[i]
            if i == n - 1:
                tot = c + madd.top_cost(S, P, k)
                if (0, 0) not in new or tot < new[(0, 0)][0]:
                    new[(0, 0)] = (tot, hist + [("TOP",)])
                continue
            for (es, ep), (cc, ops) in full_options(S, P, k).items():
                tot = c + cc
                key = (es, ep)
                if key not in new or tot < new[key][0]:
                    new[key] = (tot, hist + [ops])
        states = new
    return states[(0, 0)]


def build1(x, rows, K=0, pairs=(), cache=None):
    """single sum through madd.build_sum with a weighted plan."""
    cache = cache if cache is not None else {}
    cols, pcols, K2 = prep(rows, K & M32, N, pairs)
    key = (tuple(len(c) for c in cols), tuple((K2 >> i) & 1 for i in range(N)), tuple(len(c) for c in pcols))
    if key not in cache:
        cache[key] = plan1(list(key[0]), list(key[1]), N, list(key[2]))
    return madd.build_sum(x, rows, K & M32, plan_cache=cache, pairs=pairs)


def apply_ops(c, ops, S, P, k, es, ep):
    """apply cells (same code as madd.build_sum) to the item lists S (singles), P (pairs); returns k."""
    for op in ops:
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
        elif op in madd.EXTRA:
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
    return k


def _xor_all(c, k, S, P, acc=None):
    """XOR of constant k, singles and pair second components (top column)."""
    acc = k if acc is None else (acc ^ 1 if k else acc)
    for l in list(S) + [w for (_, w) in P]:
        if acc in (0, 1): acc = l ^ acc
        elif l in (0, 1): acc = acc ^ l
        else: acc = c.X(acc, l)
    return acc


def groups(shared, spairs, arows, erows, KS=0, KA=0, KE=0):
    gS = prep(shared, KS & M32, N, spairs)
    gA = prep(arows, KA & M32, N, ())
    gE = prep(erows, KE & M32, N, ())
    raws = tuple((len(gS[0][i]), len(gS[1][i]), len(gA[0][i]), len(gA[1][i]), len(gE[0][i]), len(gE[1][i]))
                 for i in range(N))
    return gS, gA, gE, raws, (gS[2] + gA[2]) & M32, (gS[2] + gE[2]) & M32


def estimate2(shared, spairs, arows, erows, KS=0, KA=0, KE=0, cache=None, joint=None):
    gS, gA, gE, raws, kA, kE = groups(shared, spairs, arows, erows, KS, KA, KE)
    key = (raws, kA, kE, joint)
    if cache is not None and key in cache: return cache[key][0]
    r = plan2(raws, kA, kE, joint=joint)
    if cache is not None: cache[key] = r
    return r[0]


def build2(x, shared, spairs, arows, erows, KS=0, KA=0, KE=0, cache=None, joint=None, info=None):
    """returns (A, E) = (sum(shared)+sum(arows)+KS+KA, sum(shared)+sum(erows)+KS+KE) mod 2^32."""
    c = Cells(x)
    gS, gA, gE, raws, kAtot, kEtot = groups(shared, spairs, arows, erows, KS, KA, KE)
    key = (raws, kAtot, kEtot, joint)
    if cache is not None and key in cache:
        cost, hist = cache[key]
    else:
        cost, hist = plan2(raws, kAtot, kEtot, joint=joint)
        if cache is not None: cache[key] = (cost, hist)
    if info is not None:
        info["cost"] = cost
        info["cS"] = sum(h[0] << i for i, h in enumerate(hist))
        info["rem"] = hist      # filled with remainder sizes below
    outA, outE = [], []
    cS = ([], []); cA = ([], []); cE = ([], [])
    for i in range(N):
        kS, kA, kE = hist[i][0:3]
        S = list(gS[0][i]) + cS[0]; P = list(gS[1][i]) + cS[1]
        if i == N - 1:
            sh = _xor_all(c, kS, S, P) if (S or P) else kS
            # A: shared part first, E: own part first (same operand order as column-compression generator/constant-redistribution generator)
            outA.append(_xor_all(c, kA, list(gA[0][i]) + cA[0], list(gA[1][i]) + cA[1], acc=sh))
            accE = _xor_all(c, kE, list(gE[0][i]) + cE[0], list(gE[1][i]) + cE[1])
            outE.append(_xor_all(c, 0, [sh], [], acc=accE))
            break
        opsS, kpass, opsA, opsE = hist[i][3:7]
        nS = ([], [])
        k = apply_ops(c, opsS, S, P, kS, nS[0], nS[1])
        assert k == kpass, (i, k, kpass)
        if info is not None: info.setdefault("r", []).append(len(S) + 2 * len(P) + k)
        nA = ([], []); nE = ([], [])
        SA = list(S) + list(gA[0][i]) + cA[0]; PA = list(P) + list(gA[1][i]) + cA[1]
        k2 = apply_ops(c, opsA, SA, PA, kA or kpass, nA[0], nA[1])
        assert not PA and len(SA) + k2 <= 1, (i, SA, PA, k2)
        outA.append(SA[0] if SA else k2)
        SE = list(gE[0][i]) + list(S) + cE[0]; PE = list(gE[1][i]) + list(P) + cE[1]
        k2 = apply_ops(c, opsE, SE, PE, kE or kpass, nE[0], nE[1])
        assert not PE and len(SE) + k2 <= 1, (i, SE, PE, k2)
        outE.append(SE[0] if SE else k2)
        cS, cA, cE = nS, nA, nE
    if info is not None: info["rem"] = info.pop("r", [])
    return outA, outE
