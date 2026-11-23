"""negctl.py: negative controls for the lemma pipeline.  Mutate vnew (flip an operand complement of a twinned gate, or
the lemma polarity) and check that the owning lemma is NOT proven equivalent.  usage: negctl.py COUNT SEED"""
import common, os, sys, json, random, copy
import lemlib3 as L
import equiv
from xag import XAG
REPO = os.path.join(common.HERE, "..", "..")
RP = os.environ["CERT_REF"]; NP = os.environ["CERT_NEW"]   # circuit paths from the environment (set by run_certificate.sh)
cnt, seed = int(sys.argv[1]), int(sys.argv[2])
rnd = random.Random(seed)
R = L.load(RP); N = L.load(NP)
P = json.load(open(common.HERE + "/work/plan.json"))
T = {int(k): tuple(v) for k, v in P["T"].items()}
items = [tuple(i) for i in P["items"]]
nb = N.n_inputs + 1
stat = {"not_equivalent": 0, "other": 0}
for t in range(cnt):
    u, tt, pol, d = rnd.choice(items)
    kind = rnd.choice(["pol", "flipA", "flipB"])
    N2 = XAG(N.n_inputs); N2.nodes = list(N.nodes); N2.outputs = list(N.outputs)
    it = (u, tt, pol, d)
    if kind == "pol":
        it = (u, tt, pol ^ 1, d)
    else:
        i = u - nb; op, a, b = N2.nodes[i]
        N2.nodes[i] = (op, a ^ 1, b) if kind == "flipA" else (op, a, b ^ 1)
    A, B, free = L.build(R, N2, T, [it])
    m = equiv.check(A, B, f"{common.HERE}/work/neg/n{seed}_{t}", name="neg", timeout=60)
    ok = m["status"] == "not_equivalent"
    stat["not_equivalent" if ok else "other"] += 1
    print(t, kind, "AND" if N.nodes[u - nb][0] == "AND" else "XOR", m["status"], flush=True)
print(stat)
