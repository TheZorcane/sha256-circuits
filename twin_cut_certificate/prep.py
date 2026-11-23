"""prep.py: search for plane matching + twin pairs + per-item expansion depth; writes work/plan.json
(untrusted search result: every claim in it is re-proved by SAT / re-derived by compose.py)."""
import common, json, os, sys, time, collections
import lemlib3 as L
REPO = os.path.join(common.HERE, "..", "..")
RP = os.environ["CERT_REF"]; NP = os.environ["CERT_NEW"]   # circuit paths from the environment (set by run_certificate.sh)
DMAX = 12


def pick(R, N, T, u, t, pol):
    for d in range(1, DMAX + 1):
        A, B, free = L.build(R, N, T, [(u, t, pol, d)])
        if not L.sim_equal(A, B, rounds=2):
            return d, len(A.nodes) + len(B.nodes)
    return None, None


def main():
    R = L.load(RP); N = L.load(NP)
    t0 = time.time(); M, fr, fn = L.match(R, N); print("match", len(M), round(time.time() - t0, 1), flush=True)
    T, outs = L.twins(R, N, fr, fn)
    nb = N.n_inputs + 1
    print("vnew nodes", len(fn) - nb, "twinned", len(T) - N.n_inputs, "untwinned", len(fn) - nb - (len(T) - N.n_inputs))
    items = []; dh = collections.Counter(); fail = []; size = 0
    todo = [(u, T[u][0], T[u][1]) for u in sorted(T) if u >= nb]
    todo_out = [(u, t, p) for (o, u, t, p) in outs]
    for k, (u, t, p) in enumerate(todo):
        d, sz = pick(R, N, T, u, t, p)
        if d is None:
            fail.append((u, t, p)); continue
        dh[d] += 1; size = max(size, sz)
        items.append((u, t, p, d))
        if k % 10000 == 0:
            print(k, dict(dh), round(time.time() - t0, 1), flush=True)
    oitems = []
    for (u, t, p) in todo_out:
        d, sz = pick(R, N, T, u, t, p)
        if d is None:
            fail.append((u, t, p)); continue
        dh[("o", d)] += 1
        oitems.append((u, t, p, d))
    print("depth histogram", dict(dh), "max pair size", size, "FAILED", len(fail), fail[:10])
    json.dump({"M": {str(k): v for k, v in M.items()}, "T": {str(k): v for k, v in T.items()},
               "items": items, "oitems": oitems, "outs": outs, "fail": fail}, open(common.HERE + "/work/plan.json", "w"))


main()
