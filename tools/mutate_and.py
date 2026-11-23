"""Negative control for the compositional prover: complement one input of one AND gate.
usage: mutate_and.py CAND.xag OUT.xag [round]
Picks the lowest-numbered AND gate in the interior of the cone of lemma round:R (default 10) of CAND,
complements its first operand, and writes the mutant.  lemmas.py must then report not_equivalent for that
lemma.  Run from the artifact root with PYTHONPATH=verifier."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "compositional_proof"))
import lemmas  # noqa: E402
from xag import AND, write_xag  # noqa: E402

cand, out = sys.argv[1], sys.argv[2]
rnd = int(sys.argv[3]) if len(sys.argv) > 3 else 10
x, cm, nodeslots = lemmas.get_ctx(cand)
_, _, info = lemmas.build_lemma(x, cm, nodeslots, ("round", rnd))
n_in = x.n_inputs
node = min(n for n in info["interior"] if x.nodes[n - n_in - 1][0] == AND)
op, a, b = x.nodes[node - n_in - 1]
x.nodes[node - n_in - 1] = (op, a ^ 1, b)
open(out, "w").write(write_xag(x))
print(f"mutated AND node {node} (round:{rnd} cone): first operand {a} -> {a ^ 1}")
