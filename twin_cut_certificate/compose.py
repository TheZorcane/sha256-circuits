#!/usr/bin/env python3
"""compose.py WORKDIR : composition checker for the twin-cut certificate.

Re-derives everything needed for the induction from the raw .xag files and the batch records, with its OWN lemma
builder (recursive/memoised, different from lemlib3.build), and checks:
 C1  every batch manifest: status equivalent, both LRAT checkers accepted, proof file hash == manifest, CNF text
     recomputed from the rebuilt circuits hashes to the manifest's cnf_sha256 (the proof is about THESE circuits)
 C2  rebuilt circuits (lemlib3, the producer) have the manifest circuit hashes; independent builder's circuits are
     structurally identical to them (canonical form hash equal)
 C3  twin table T (from the batch items; inputs are twins of themselves) is a function; every vnew A-leaf s of every
     lemma is an input or has a lemma (s,T[s]) in the set and s < u  (=> induction in vnew order is closed)
 C4  all 256 output lemmas present in SPEC order: for output o, u = N.outputs[o]>>1, t = R.outputs[o]>>1 and
     (N.outputs[o]&1) ^ pol == (R.outputs[o]&1); the lemma for (u,t,pol) is in the proven set
 C5  coverage statistics (vnew nodes twinned, interior non-twinned, dead nodes)
 C6  semantic faithfulness by simulation (2048 random vectors): value_new(u)^pol == value_ref(t) for all items, and each
     batch's circuits A/B evaluated at value_ref(free vars) give value_new(u)^pol / value_ref(t)
Induction: for every vnew node s with lemma: value_new(s) == value_ref(T[s]) ^ pol (by lemma at the true values of the
free variables, A-leaves true by induction in vnew order, vref vars are vref node values by construction).  For outputs the same
for the output lemma, hence N(x) == R(x).
"""
import common, sys, os, json, random, hashlib, glob
sys.setrecursionlimit(100000)
import equiv, cost
from xag import XAG, AND, XOR, read_xag
import lemlib3 as L   # ONLY used as the producer whose raw circuits the manifests were proven for

REPO = os.path.join(common.HERE, "..", "..")
RP = os.environ["CERT_REF"]; NP = os.environ["CERT_NEW"]   # circuit paths from the environment (set by run_certificate.sh)


def indep_build(R, N, T, items):
    rb, nb = R.n_inputs + 1, N.n_inputs + 1
    def r_ops(v): return R.nodes[v - rb]
    info = []
    for (u, t, pol, d) in items:
        # A-leaves
        leaves8 = set(); seen = set()
        def walk(v, top):
            if v == 0 or v in seen and not top:
                return
            if not top and v in T:
                leaves8.add(v); return
            seen.add(v)
            assert v >= nb
            _, a, b = N.nodes[v - nb]
            walk(a >> 1, False); walk(b >> 1, False)
        walk(u, True)
        roots = {T[s][0] for s in leaves8} | {t}
        # X = vref gates within distance < d of a root
        level = {r for r in roots if r >= rb}; X = set()
        for _ in range(d):
            X |= level
            level = {x >> 1 for v in level for x in r_ops(v)[1:] if (x >> 1) >= rb}
        refs = set()
        for r in roots: refs.add(r)
        for v in X:
            for x in r_ops(v)[1:]: refs.add(x >> 1)
        fv = {v for v in refs if v and (v < rb or v not in X)}
        info.append((u, t, pol, leaves8, X, fv))
    free = sorted(set().union(*[i[5] for i in info]))
    pos = {v: i for i, v in enumerate(free)}
    A = XAG(len(free)); B = XAG(len(free))
    def mk(C):
        memo = {}
        def g7(v, X):
            if v in memo: return memo[v]
            if v not in X:
                r = 2 * (pos[v] + 1)
            else:
                op, a, b = r_ops(v)
                la, lb = g7(a >> 1, X) ^ (a & 1), g7(b >> 1, X) ^ (b & 1)
                r = C.AND(la, lb) if op == AND else C.XOR(la, lb)
            memo[v] = r; return r
        return memo, g7
    oa, ob = [], []
    for (u, t, pol, leaves8, X, fv) in info:
        memoA, gA = mk(A); memoB, gB = mk(B)
        m8 = {}
        def g8(v, top):
            if not top and v in T:
                tt, p = T[v]; return gA(tt, X) ^ p
            if v in m8: return m8[v]
            op, a, b = N.nodes[v - nb]
            la, lb = g8(a >> 1, False) ^ (a & 1), g8(b >> 1, False) ^ (b & 1)
            r = A.AND(la, lb) if op == AND else A.XOR(la, lb)
            m8[v] = r; return r
        oa.append(g8(u, True) ^ pol); ob.append(gB(t, X))
    A.outputs = oa; B.outputs = ob
    return A, B, free


def main(work):
    R = read_xag(open(RP).read()); N = read_xag(open(NP).read())
    rb, nb = R.n_inputs + 1, N.n_inputs + 1
    n = R.n_inputs
    assert R.n_inputs == N.n_inputs == 512 and len(R.outputs) == len(N.outputs) == 256
    err = []; E = err.append
    batches = []
    for d in sorted(glob.glob(work + "/b*/batch.json")):
        batches.append(json.load(open(d)))
    allitems = []
    for b in batches:
        allitems += [tuple(i) for i in b["items"]]
    # which batch holds outputs? outputs are items flagged by being the final batches: recompute by identity below
    P = json.load(open(common.HERE + "/work/plan.json"))
    out_items = [tuple(i) for i in P["oitems"]]
    gate_items = [tuple(i) for i in P["items"]]
    assert sorted(allitems) == sorted(gate_items + out_items), "batch items differ from plan"
    # twin table from gate items only (inputs identity)
    T = {i: (i, 0) for i in range(1, nb)}
    for (u, t, pol, d) in gate_items:
        if u in T:
            E(f"C3 duplicate lemma for node {u}")
        T[u] = (t, pol)
    # C1/C2
    cnt_free = 0; cnt_cl = 0; proof_b = 0
    for b in batches:
        name = b["name"]; wd = f"{work}/{name}"
        m = json.load(open(f"{wd}/{name}.manifest.json"))
        pc = m["proof_check"]
        if not (m["status"] == "equivalent" and pc["accepted"] and pc["python_checker"]["accepted"] and pc["vendored_lrat_check"]["accepted"]):
            E(f"C1 {name} proof not accepted")
        if os.path.exists(f"{wd}/{name}.lrat"):
            if equiv.sha256_file(f"{wd}/{name}.lrat") != m["proof_sha256"]:
                E(f"C1 {name} proof file hash mismatch")
        items = [tuple(i) for i in b["items"]]
        A, B, free = L.build(R, N, T, items)
        if free != b["free"]: E(f"C2 {name} free list differs")
        if cost.xag_hash(A) != m["circuit_a_sha256"] or cost.xag_hash(B) != m["circuit_b_sha256"]:
            E(f"C2 {name} raw circuit hash differs from proven manifest")
        cnf_text, _ = equiv.miter_cnf(A, B)
        if hashlib.sha256(cnf_text.encode()).hexdigest() != m["cnf_sha256"]:
            E(f"C1 {name} CNF hash differs from proven CNF")
        A2, B2, free2 = indep_build(R, N, T, items)
        if free2 != free or cost.xag_hash(cost.canonicalize(A2)) != cost.xag_hash(cost.canonicalize(A)) \
                or cost.xag_hash(cost.canonicalize(B2)) != cost.xag_hash(cost.canonicalize(B)):
            E(f"C2 {name} independent rebuild differs structurally")
        cnt_free += len(free); proof_b += m.get("proof_bytes") or 0
        b["_AB"] = (A, B)
    # C3 leaf closure (independent walk)
    def a_leaves(u):
        out = set(); seen = set()
        def walk(v, top):
            if v == 0: return
            if not top and v in T:
                out.add(v); return
            if v in seen: return
            seen.add(v); assert v >= nb
            _, a, b = N.nodes[v - nb]
            walk(a >> 1, False); walk(b >> 1, False)
        walk(u, True)
        return out
    for (u, t, pol, d) in gate_items + out_items:
        for s in a_leaves(u):
            if not (s < nb or s in T):
                E(f"C3 leaf {s} of lemma {u} has no lemma")
            if s >= u:
                E(f"C3 leaf {s} >= {u}")
    # C4 outputs
    for o in range(256):
        lo, ro = N.outputs[o], R.outputs[o]
        u, t = lo >> 1, ro >> 1
        its = [i for i in out_items if i[0] == u and i[1] == t]
        if u < nb or u == 0 or t == 0:
            E(f"C4 output {o} node is input/const (unexpected)"); continue
        if not its:
            E(f"C4 output {o} lemma missing"); continue
        pol = its[0][2]
        if ((lo & 1) ^ pol) != (ro & 1):
            E(f"C4 output {o} polarity mismatch")
    # C5 coverage
    reach = set(); stack = [o >> 1 for o in N.outputs]
    while stack:
        v = stack.pop()
        if v < nb or v in reach: continue
        reach.add(v); _, a, b = N.nodes[v - nb]; stack += [a >> 1, b >> 1]
    dead = (len(N.nodes)) - len(reach)
    twinned = sum(1 for v in range(nb, nb + len(N.nodes)) if v in T)
    # C6 simulation
    W = 2048; mask = (1 << W) - 1; rnd = random.Random(0xD06)
    lns_ = [rnd.getrandbits(W) for _ in range(n)]
    def allvals(X):
        val = [0] + lns_[:]
        for op, a, b in X.nodes:
            x = val[a >> 1] ^ (mask if a & 1 else 0); y = val[b >> 1] ^ (mask if b & 1 else 0)
            val.append(x & y if op == AND else x ^ y)
        return val
    vref, vnew = allvals(R), allvals(N)
    for (u, t, pol, d) in gate_items + out_items:
        if vnew[u] ^ (mask if pol else 0) != vref[t]:
            E(f"C6 twin ({u},{t},{pol}) false in simulation")
    for b in batches:
        A, B = b["_AB"]; vals = [0] + [vref[v] for v in b["free"]]
        for X, which in ((A, "A"), (B, "B")):
            val = list(vals)
            for op, a, bb in X.nodes:
                x = val[a >> 1] ^ (mask if a & 1 else 0); y = val[bb >> 1] ^ (mask if bb & 1 else 0)
                val.append(x & y if op == AND else x ^ y)
            outs = [val[o >> 1] ^ (mask if o & 1 else 0) for o in X.outputs]
            for (u, t, pol, d), ov in zip(b["items"], outs):
                want = (vnew[u] ^ (mask if pol else 0)) if which == "A" else vref[t]
                if ov != want:
                    E(f"C6 batch {b['name']} {which} item {u} not faithful in simulation"); break
    rep = {"batches": len(batches), "gate_lemmas": len(gate_items), "output_lemmas": len(out_items),
           "new_gate_nodes": len(N.nodes), "twinned_gate_nodes": twinned, "interior_untwinned": len(N.nodes) - twinned,
           "dead_nodes": dead, "total_free_vars_over_batches": cnt_free, "proof_bytes_total": proof_b, "errors": err[:20], "n_errors": len(err)}
    print(json.dumps(rep, indent=1))
    return 0 if not err else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
