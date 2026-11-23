"""Step 3: per-lemma cone extraction + proof against the verifier's reference (sharef.py).
usage: lemmas.py CAND.xag OUTDIR [--jobs 5] [--only kind:i,...] [--timeout 1800]
Lemma = (kind,i) in sharef.lemma_list().  For each lemma:
  * input slots (name,j) of the declared input words; constants (lits 0/1) stay constants,
    other slots become the free inputs of the miter.  The cut map (cutmap.build_cutmap) gives ONE signed
    literal per slot (the most-used signature match; other matches are ordinary interior nodes).
  * a polarity-aware union-find over candidate nodes gives one free variable per distinct leaf node; a slot's
    input literal is that variable XOR the slot literal's sign, so two slots served by the same node share a
    variable (with the right polarity).  The code would also merge several literals of one slot, but the cut
    map never supplies more than one, and check_compose.py (C2) rejects an output slot produced twice.
  * round 4 (EXPAND): a_0 and e_0 are not leaves; the cone is traversed through the candidate's own logic down
    to primary inputs (the bits of W_0), and the reference receives the same candidate expansion of a_0, e_0.
  * cone A = sub-circuit from the leaves to the literals of the output slots.
  * reference B = sharef.lemma_apply (style A) on the same signed input words, in bit order.
  * a miter equiv.check A vs B must be `equivalent`.
"""
import common, sys, json, os, time, argparse
from multiprocessing import Pool
import cost, equiv
from xag import XAG, AND, XOR, read_xag
import sharef
from sharef import lemma_list, lemma_in_out, lemma_apply, iv_word


class ConeError(Exception):
    pass


class MissingSlot(ConeError):
    """No literal of the candidate matches a reference word bit on the cut-map simulation vectors."""
    pass


def slot_consts(name):
    v = iv_word(name)
    return None if v is None else sharef.cword(v)


def prepare(x, cm):
    """node -> list of (slotname, j) over all slots (for foreign-cut detection)"""
    nodeslots = {}
    # column-compression generator (copy only): with VERIFY_TAIL=1 only slots that are some lemma's input/output count as cuts
    used = None
    if sharef.TAIL:
        used = set()
        for l in lemma_list():
            i_, o_ = lemma_in_out(l)
            used |= set(i_) | set(o_)
    for (name, j), alts in cm.items():
        if used is not None and name not in used:
            continue
        for L in alts:
            nodeslots.setdefault(L >> 1, []).append((name, j))
    return nodeslots


EXPAND = {('round', 4): {('a', 0), ('e', 0)}}


def build_lemma(x, cm, nodeslots, l, expand=None):
    ins, outs = lemma_in_out(l)
    expand = EXPAND.get(l, set()) if expand is None else expand
    expnames = [n for n in ins if n in expand]
    # input slots
    # weighted union-find over candidate nodes: value(node) = value(root) ^ parity.  Nodes are merged
    # when they are alternatives of one slot, or the same node serves two slots (aliasing).
    par = {}

    def find(n):
        if n not in par:
            par[n] = (n, 0)
        r, p = par[n]
        if r == n:
            return n, 0
        rr, pp = find(r)
        par[n] = (rr, p ^ pp)
        return rr, p ^ pp

    def union(n1, n2, rel):   # value(n1) ^ value(n2) == rel
        r1, p1 = find(n1); r2, p2 = find(n2)
        if r1 == r2:
            if p1 ^ p2 != rel:
                raise ConeError(f"inconsistent aliasing between nodes {n1} and {n2}")
            return
        par[r1] = (r2, p1 ^ p2 ^ rel)
    slotnodes = {}
    inw = {}
    for name in ins:
        cw = slot_consts(name)
        if cw is not None:
            inw[name] = cw
            continue
        if name in expand:
            continue
        for j in range(32):
            alts = cm.get((name, j)) or []
            if not alts:
                raise MissingSlot(f"input word bit {name}[{j}] not found in the candidate")
            n0, neg0 = alts[0] >> 1, alts[0] & 1
            slotnodes[(name, j)] = (n0, neg0)
            for L in alts[1:]:
                union(L >> 1, n0, (L & 1) ^ neg0)
            for L in alts:
                find(L >> 1)
    roots = sorted({find(n)[0] for n in list(par)})
    ridx = {r: i for i, r in enumerate(roots)}
    slots = [None] * len(roots)
    leaf = {}
    for n in list(par):
        r, p = find(n)
        leaf[n] = (ridx[r], p)
    for name in ins:
        if name in inw or name in expand:
            continue
        w = []
        for j in range(32):
            n0, neg0 = slotnodes[(name, j)]
            r, p = find(n0)
            w.append(2 * (ridx[r] + 1) ^ p ^ neg0)
            slots[ridx[r]] = slots[ridx[r]] or (name, j)
        inw[name] = w
    ns = len(slots)
    # outputs: all alternatives
    outlits, outslot = [], []
    for name in outs:
        for j in range(32):
            if not cm.get((name, j)):
                raise MissingSlot(f"output word bit {name}[{j}] not found in the candidate")
            for L in cm[(name, j)]:
                outlits.append(L); outslot.append((name, j))
    own = {L >> 1 for L in outlits}
    # collect interior
    n_in = x.n_inputs
    pis = set()

    def trav(starts):
        inter, stack = set(), list(starts)
        while stack:
            n = stack.pop()
            if n == 0 or n in leaf or n in inter:
                continue
            if n <= n_in:
                if not expand:
                    raise ConeError(f"hidden dependency on primary input {n - 1}")
                pis.add(n); continue
            if n in nodeslots and n not in own and n not in expnodes:
                raise ConeError(f"hidden dependency on foreign cut node {n} {nodeslots[n][:2]}")
            inter.add(n)
            op, a, b = x.nodes[n - n_in - 1]
            stack.append(a >> 1); stack.append(b >> 1)
        return inter
    expnodes = {cm[(nm, j)][0] >> 1 for nm in expnames for j in range(32)}
    inter = trav([L >> 1 for L in outlits])
    inter2 = trav(list(expnodes)) if expnames else set()
    pis = sorted(pis)
    nt = ns + len(pis)
    mp = {0: 0}
    for n, (idx, neg) in leaf.items():
        mp[n] = 2 * (idx + 1) ^ neg
    for k, n in enumerate(pis):
        mp[n] = 2 * (ns + k + 1)

    def build_into(xg, nodes):
        def m(l): return mp[l >> 1] ^ (l & 1)
        for n in sorted(nodes):
            op, a, b = x.nodes[n - n_in - 1]
            mp[n] = xg.AND(m(a), m(b)) if op == AND else xg.XOR(m(a), m(b))
        return m
    cone = XAG(nt)
    m = build_into(cone, inter)
    cone.outputs = [m(L) for L in outlits]
    # reference (expanded slots re-use the candidate's own expansion down to primary inputs; sound because
    # the producing lemma proves those nodes equal the reference words as functions of the primary inputs)
    ref = XAG(nt)
    for n in list(mp):
        if n in inter or n in inter2:
            del mp[n]
    mr = build_into(ref, inter | inter2)
    for nm in expnames:
        inw[nm] = [mr(cm[(nm, j)][0]) for j in range(32)]
    res = lemma_apply(ref, l, inw, "A")
    byname = dict(zip(outs, res))
    ref.outputs = [byname[nm][j] for (nm, j) in outslot]
    ns = nt
    info = {"lemma": f"{l[0]}:{l[1]}", "n_free": ns, "expanded": [list(n) for n in expnames], "pi_free": [p - 1 for p in pis], "n_interior": len(inter),
            "interior": sorted(inter), "n_out_lits": len(outlits),
            "out_lits": outlits, "out_slots": [(list(a), b) for a, b in outslot],
            "in_names": [list(n) for n in ins], "leaf_nodes": sorted(leaf)}
    return cone, ref, info


_CTX = {}


def get_ctx(candpath):
    if candpath not in _CTX:
        import pickle
        from cutmap import build_cutmap
        x = read_xag(open(candpath).read())
        pk = candpath + ".cm.pkl"
        if os.path.exists(pk):
            cm = pickle.load(open(pk, "rb"))
        else:
            cm = build_cutmap(x, None)
            pickle.dump(cm, open(pk, "wb"))
        _CTX[candpath] = (x, cm, prepare(x, cm))
    return _CTX[candpath]


_SIM = {}


def whole_sim(candpath, n=2000):
    if candpath not in _SIM:
        import sim
        x = get_ctx(candpath)[0]
        _SIM[candpath] = sim.verify(x, "primary", list(range(512)), [(s, 0) for s in range(256)], n=n, seed=1)
    return _SIM[candpath]


def prove_one(args):
    candpath, outdir, l, timeout = args
    x, cm, nodeslots = get_ctx(candpath)
    name = f"{l[0]}_{l[1]}"
    try:
        cone, ref, info = build_lemma(x, cm, nodeslots, l)
    except MissingSlot as e:
        # A word bit of the reference is computed by no node of the candidate (on the cut-map simulation
        # vectors).  Report not_equivalent only if a whole-circuit simulation against the FIPS 180-4
        # reference finds a confirmed counterexample; otherwise this is a cone error (never a pass).
        sv = whole_sim(candpath)
        if sv["mismatches"]:
            return {"lemma": f"{l[0]}:{l[1]}", "status": "not_equivalent", "error": str(e),
                    "evidence": "whole-circuit simulation against the reference",
                    "first_mismatch": sv["first_mismatch"]}
        return {"lemma": f"{l[0]}:{l[1]}", "status": "cone_error", "error": str(e)}
    except ConeError as e:
        return {"lemma": f"{l[0]}:{l[1]}", "status": "cone_error", "error": str(e)}
    t = time.time()
    mf = equiv.check(cone, ref, os.path.join(outdir, name), name=name, timeout=timeout,
                     interface="free inputs = non-constant input word-bit slots; outputs = all alt literals of output slots")
    info.update({"status": mf["status"], "wall_s": round(time.time() - t, 2),
                 "cone_hash": mf["circuit_a_sha256"], "ref_hash": mf["circuit_b_sha256"],
                 "cnf_sha256": mf["cnf_sha256"], "proof_sha256": mf.get("proof_sha256"),
                 "cone_counts": mf["circuit_a_counts"], "ref_counts": mf["circuit_b_counts"]})
    # drop bulky proof artifacts after verification (hashes recorded in manifest)
    for ext in (".lrat", ".cnf"):
        p = os.path.join(outdir, name, name + ext)
        if info["status"] == "equivalent" and os.path.exists(p):
            os.remove(p)
    return info


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cand"); ap.add_argument("outdir")
    ap.add_argument("--jobs", type=int, default=5)
    ap.add_argument("--only", default="")
    ap.add_argument("--timeout", type=float, default=1800)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    L = lemma_list()
    if a.only:
        want = {tuple((p.split(":")[0], int(p.split(":")[1]))) for p in a.only.split(",")}
        L = [l for l in L if l in want]
    t0 = time.time()
    get_ctx(a.cand)
    with Pool(a.jobs) as pool:
        results = []
        for r in pool.imap_unordered(prove_one, [(a.cand, a.outdir, l, a.timeout) for l in L]):
            print(r["lemma"], r["status"], r.get("wall_s"), r.get("error", ""), flush=True)
            results.append(r)
    results.sort(key=lambda r: r["lemma"])
    json.dump(results, open(os.path.join(a.outdir, "lemmas.json"), "w"))
    from collections import Counter
    print(Counter(r["status"] for r in results), round(time.time() - t0, 1), "s")
