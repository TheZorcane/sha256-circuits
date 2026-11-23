"""Step 5: attacks.  (a) edge-pattern simulation of the candidate vs trusted sha256_ref;
(b) mutation test: single-gate mutants of the candidate must be caught by the lemma pipeline
(mutant keeps node numbering so the ORIGINAL cut map / interface claim is reused)."""
import common, sys, json, os, random, shutil, subprocess, time
import sim, cost
from xag import read_xag, write_xag, AND, XOR
import sha256_ref as S

x = read_xag(open("cand.xag").read())
rep = {}
# ---------------- (a) edge simulation ----------------
M = S.MASK
vs = [[0] * 16, [M] * 16]
for t in range(16):
    for b in range(32):
        w = [0] * 16; w[t] = 1 << b; vs.append(w)            # single set bit (512)
        w = [M] * 16; w[t] = M ^ (1 << b); vs.append(w)       # single cleared bit (512)
pats = [0x55555555, 0xaaaaaaaa, 0x80000000, 0x7fffffff, 0x00000001, 0xfffffffe, 0x0000ffff, 0xffff0000,
        0x01010101, 0x80808080]
for p in pats:
    vs.append([p] * 16)
for p in pats:
    for q in pats:
        vs.append([p, q] * 8)
for k in range(1, 33):                     # carries rippling to/through bit 31: words 2^k-1 and complements
    vs.append([(1 << k) - 1] * 16); vs.append([M ^ ((1 << k) - 1)] * 16)
    vs.append([(1 << k) - 1 if i % 2 else 1 for i in range(16)])
for m in (b"", b"abc", b"a" * 55, b"a" * 56 and b"a" * 55):
    vs.append(S.pad(m)[0])
r = random.Random(31337)
for _ in range(20000):                      # sparse / dense random
    d = r.choice([0.02, 0.5, 0.98])
    vs.append([sum((1 if r.random() < d else 0) << b for b in range(32)) for _ in range(16)])
N = len(vs)
lns_ = sim.vectors_to_lanes([S.words_to_int(w) for w in vs], 512)
out = sim.lanes_to_vectors(sim.simulate(x, lns_, N), N)
bad = sum(1 for o, w in zip(out, vs) if o != S.words_to_int(S.compress_fixed_iv(w)))
rep["edge_sim"] = {"vectors": N, "mismatches": bad}
print(rep["edge_sim"], flush=True)

# ---------------- (b) mutation test ----------------
n_in = x.n_inputs
R3 = {r["lemma"]: r for r in json.load(open("work/run3/lemmas.json"))}
own = {}
for k, rr in R3.items():
    for n in rr["interior"]:
        own.setdefault(n, []).append(k)
rng = random.Random(2026)
# targeted picks: one node in schedule, early round, mid round, late round, final, plus randoms
def pick(prefix):
    c = [n for n, v in own.items() if any(s.startswith(prefix) for s in v)]
    return rng.choice(sorted(c))
targets = [pick("sched:20"), pick("sched:63"), pick("round:0"), pick("round:2"), pick("round:30"), pick("round:63"),
           pick("final:0"), pick("final:7")]
# a bit-31 carry-related node: highest-index node of final:3 / round:31 interior
targets.append(max(R3["final:3"]["interior"])); targets.append(max(R3["round:31"]["interior"]) - 5)
targets += [rng.randint(n_in + 1, n_in + len(x.nodes)) for _ in range(4)]
kinds = ["op", "edge", "op", "edge", "op", "edge", "op", "edge", "op", "edge", "op", "edge"]
results = []
for ti, (n, kind) in enumerate(zip(targets, kinds)):
    y = read_xag(write_xag(x))
    k = n - n_in - 1
    op, a, b = y.nodes[k]
    if kind == "op":
        y.nodes[k] = (XOR if op == AND else AND, a, b)
    else:
        y.nodes[k] = (op, a ^ 1, b) if op == AND else (op, a, b ^ 1) if False else (op, a ^ 1, b)
    # skip if mutant is trivially identical after canonical form (XOR edge flip only complements output)
    name = f"mut{ti}"
    path = f"work/{name}.xag"
    open(path, "w").write(write_xag(y))
    shutil.copy("cand.xag.cm.pkl", path + ".cm.pkl")
    # functional change visible by simulation?
    o2 = sim.lanes_to_vectors(sim.simulate(y, lns_[:], N), N)
    sim_diff = sum(1 for o, w in zip(o2, vs) if o != S.words_to_int(S.compress_fixed_iv(w)))
    t0 = time.time()
    p = subprocess.run([sys.executable, "lemmas.py", path, f"work/{name}", "--jobs", "5", "--timeout", "600"],
                       capture_output=True, text=True, env=dict(os.environ, PYTHONPATH="../../trusted"))
    last = p.stdout.strip().splitlines()[-1]
    fails = [ln for ln in p.stdout.splitlines() if " equivalent" not in ln and ln.split()[0].split(":")[0] in ("sched", "round", "final")]
    results.append({"mutant": name, "node": n, "kind": kind, "owner_lemmas": own.get(n), "sim_mismatching_vectors": sim_diff,
                    "failed_lemmas": [f.split()[0] + ":" + f.split()[1] for f in fails], "summary": last,
                    "wall_s": round(time.time() - t0, 1)})
    print(results[-1], flush=True)
rep["mutants"] = results
json.dump(rep, open("work/attack_report.json", "w"), indent=1)
