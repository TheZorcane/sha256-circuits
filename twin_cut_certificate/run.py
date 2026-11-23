"""run.py: build batched twin lemmas from work/plan.json and prove them with trusted equiv.check.
usage: run.py OUTDIR [--batch 500] [--jobs 5] [--limit K] [--only i,j]"""
import common, json, os, sys, time, argparse
from concurrent.futures import ProcessPoolExecutor
import lemlib3 as L
import equiv, cost
REPO = os.path.join(common.HERE, "..", "..")
RP = os.environ["CERT_REF"]; NP = os.environ["CERT_NEW"]   # circuit paths from the environment (set by run_certificate.sh)
G = {}


def init():
    R = L.load(RP); N = L.load(NP)
    P = json.load(open(common.HERE + "/work/plan.json"))
    T = {int(k): tuple(v) for k, v in P["T"].items()}
    G.update(R=R, N=N, T=T, TW=set(t for t, p in T.values()), P=P)


def batches(items, size):
    return [items[i:i + size] for i in range(0, len(items), size)]


def work(args):
    bi, items, outdir = args
    R, N, T, TW = G["R"], G["N"], G["T"], G["TW"]
    name = f"b{bi:04d}"
    wd = os.path.join(outdir, name)
    mp = os.path.join(wd, name + ".manifest.json")
    A, B, free = L.build(R, N, T, [tuple(i) for i in items])
    t0 = time.time()
    m = equiv.check(A, B, wd, name=name, timeout=900, extra={"batch": bi, "n_items": len(items), "free_vars": len(free)})
    json.dump({"name": name, "items": items, "free": free, "status": m["status"], "hashA": m["circuit_a_sha256"],
               "hashB": m["circuit_b_sha256"], "cnf": m["cnf_vars_clauses"], "wall": m["total_wall_s"],
               "proof_bytes": m.get("proof_bytes")}, open(os.path.join(wd, "batch.json"), "w"))
    # drop big files once accepted (cnf/lrat are regenerable; keep hashes)
    return name, m["status"], m["total_wall_s"], m.get("proof_bytes"), m["cnf_vars_clauses"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir"); ap.add_argument("--batch", type=int, default=500)
    ap.add_argument("--jobs", type=int, default=5); ap.add_argument("--limit", type=int)
    ap.add_argument("--only")
    a = ap.parse_args()
    init()
    P = G["P"]
    items = [tuple(i) for i in P["items"]]
    oit = [tuple(i) for i in P["oitems"]]
    bs = batches(items, a.batch) + batches(oit, 256)
    jobs = [(i, b, a.outdir) for i, b in enumerate(bs)]
    if a.only:
        s = set(int(x) for x in a.only.split(","))
        jobs = [j for j in jobs if j[0] in s]
    if a.limit:
        jobs = jobs[:a.limit]
    os.makedirs(a.outdir, exist_ok=True)
    t0 = time.time(); bad = 0
    with ProcessPoolExecutor(a.jobs, initializer=init) as ex:
        for name, st, w, pb, cv in ex.map(work, jobs):
            print(name, st, w, pb, cv, flush=True)
            if st != "equivalent":
                bad += 1
    print("done", len(jobs), "batches, bad", bad, "wall", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
