#!/usr/bin/env python3
"""Negative controls for the two GF(2) checkers: functionally different mutants of a final circuit must be rejected.

usage: python3 mutants.py REF.xag NEW.xag COUNT SEED
Generates COUNT mutants of NEW (operand complement flipped, AND turned into XOR or back, operand rewired to an
earlier signal, output complemented), keeps only mutants whose outputs differ from NEW on 1,024 random vectors
(so each is functionally different), and runs lincheck.py and independent/check2.py on (REF, mutant).  A checker
accepts when its last output line starts with PASS.  Prints one line per mutant and
"MUTANTS REJECTED: k of k by both checkers" when no checker accepted any mutant."""
import os
import random
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "verifier"))
import sim                               # noqa: E402
from xag import XAG, AND, XOR, read_xag, write_xag  # noqa: E402


def accepts(cmd, cwd):
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    return bool(lines) and lines[-1].startswith("PASS")


def main(refp, newp, count, seed):
    rnd = random.Random(seed)
    N = read_xag(open(newp).read())
    n, nb, W = N.n_inputs, N.n_inputs + 1, 1024
    inbits = [rnd.getrandbits(W) for _ in range(n)]

    def outs(X):
        return sim.simulate(X, inbits, W)
    base = outs(N)
    ands = [i for i, (op, a, b) in enumerate(N.nodes) if op == AND]
    xors = [i for i, (op, a, b) in enumerate(N.nodes) if op == XOR]
    tmp = tempfile.mkdtemp()
    refa = os.path.abspath(refp)
    done = rejected = 0
    while done < count:
        nodes, outputs = list(N.nodes), list(N.outputs)
        kind = rnd.choice(["flip", "and2xor", "xor2and", "rewire", "outflip"])
        if kind == "flip":
            i = rnd.randrange(len(nodes)); op, a, b = nodes[i]
            nodes[i] = (op, a ^ 1, b) if rnd.random() < .5 else (op, a, b ^ 1)
        elif kind == "and2xor":
            i = rnd.choice(ands); _, a, b = nodes[i]; nodes[i] = (XOR, a, b)
        elif kind == "xor2and":
            i = rnd.choice(xors); _, a, b = nodes[i]; nodes[i] = (AND, a, b)
        elif kind == "rewire":
            i = rnd.randrange(len(nodes)); op, a, b = nodes[i]
            c = 2 * rnd.randrange(1, nb + i) + rnd.randrange(2)
            nodes[i] = (op, c, b) if rnd.random() < .5 else (op, a, c)
        else:
            outputs[rnd.randrange(len(outputs))] ^= 1
        X = XAG(n); X.nodes = nodes; X.outputs = outputs
        try:
            X.validate()
        except Exception:
            continue
        if outs(X) == base:
            continue                     # not functionally different on the sample: skip
        p = os.path.join(tmp, "mutant_%d.xag" % done)
        open(p, "w").write(write_xag(X))
        a1 = accepts([sys.executable, os.path.join(HERE, "lincheck.py"), refa, p], HERE)
        a2 = accepts([sys.executable, "check2.py", refa, p], os.path.join(HERE, "independent"))
        ok = not a1 and not a2
        rejected += ok
        print("mutant %2d %-8s lincheck %s  check2 %s  %s" % (done, kind, "accepted" if a1 else "rejected",
              "accepted" if a2 else "rejected", "OK" if ok else "SOUNDNESS FAILURE"), flush=True)
        os.remove(p)
        done += 1
    os.rmdir(tmp)
    print(("MUTANTS REJECTED: %d of %d by both checkers" % (rejected, count)) if rejected == count
          else "MUTANT ACCEPTED: %d of %d rejected" % (rejected, count))
    return rejected == count


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    sys.exit(0 if main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])) else 1)
