"""joint-sum planner: one local rewrite, found by diffing the step-7 circuit of the paper's lineage (Table 3) against its ABC `map -a` result (the 32 XORs ABC saves):

    n = XOR(u, m),  m = AND(u^alpha, q^beta),  q = XOR(u, r),  m used only by n
    =>  n = AND(u^s1, r^s2)^s3            (u ^ (u & (u^r)) = u & r, and the three complemented variants)

The AND m is replaced by another AND (AND count unchanged), the XOR n disappears; q dies if it had no other
user.  The pattern is produced by carry cells whose pair partner degenerates in the constant-heavy first
rounds.  The rewrite is local and exact (the 4-row truth table is evaluated, not assumed); the result is
re-canonicalised.  Equivalence of the final circuit is established separately (verify/), not by this file.

usage: peep.py IN.xag OUT.xag
"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
from xag import XAG, read_xag, write_xag
import cost


def peep(x):
    base = x.n_inputs + 1
    fan = {}
    for op, a, b in x.nodes:
        for l in (a, b): fan[l >> 1] = fan.get(l >> 1, 0) + 1
    for o in x.outputs: fan[o >> 1] = fan.get(o >> 1, 0) + 1
    y = XAG(x.n_inputs)
    mp = list(range(0, 2 * base, 2))            # old var -> new literal
    def tr(l): return mp[l >> 1] ^ (l & 1)
    hits = 0
    for k, (op, a, b) in enumerate(x.nodes):
        new = None
        if op == "XOR":
            for ul, ml in ((a, b), (b, a)):
                u, m = ul >> 1, ml >> 1
                if m < base or x.nodes[m - base][0] != "AND" or fan.get(m, 0) != 1: continue
                _, ma, mb = x.nodes[m - base]
                for p, ql in ((ma, mb), (mb, ma)):
                    q = ql >> 1
                    if p >> 1 != u or q < base or x.nodes[q - base][0] != "XOR": continue
                    _, qa, qb = x.nodes[q - base]
                    if qa >> 1 == u: rl = qb
                    elif qb >> 1 == u: rl = qa
                    else: continue
                    r = rl >> 1
                    if r == u: continue
                    # truth table of n over (u, r) with all edge polarities
                    tt = []
                    for uv in (0, 1):
                        for rv in (0, 1):
                            qv = (uv ^ (1 if (qa >> 1 == u and qa & 1) or (qb >> 1 == u and qb & 1) else 0)) ^ (rv ^ (rl & 1))
                            mv = (uv ^ (p & 1)) & (qv ^ (ql & 1))
                            tt.append((uv ^ (ul & 1)) ^ (mv ^ (ml & 1)))
                    ones = sum(tt)
                    if ones not in (1, 3): continue
                    s3 = 1 if ones == 3 else 0
                    idx = [i for i, v in enumerate(tt) if v != s3][0]      # the single minterm of n^s3
                    s1, s2 = 1 - (idx >> 1), 1 - (idx & 1)
                    new = y.AND(mp[u] ^ s1, mp[r] ^ s2) ^ s3
                    hits += 1
                    break
                if new is not None: break
        if new is None:
            new = y.AND(tr(a), tr(b)) if op == "AND" else y.XOR(tr(a), tr(b))
        mp.append(new)
    y.outputs = [tr(o) for o in x.outputs]
    return cost.canonicalize(y), hits


if __name__ == "__main__":
    x = read_xag(open(sys.argv[1]).read())
    tot = 0
    while True:
        y, h = peep(x)
        tot += h
        if h == 0 or y.counts()["total"] >= x.counts()["total"]:
            if y.counts()["total"] < x.counts()["total"]: x = y
            break
        x = y
    r = cost.report_xag(x); r["peephole_hits"] = tot
    print(json.dumps(r))
    open(sys.argv[2], "w").write(write_xag(x))
