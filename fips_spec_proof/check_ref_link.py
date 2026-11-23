#!/usr/bin/env python3
"""Check that the reference side of every SAT lemma is node for node the fragment proved in SAW.

usage: python3 check_ref_link.py CIRCUIT.xag VARIANT [BUILD_DIR]
       VARIANT = split (the two proof bases) | tail (the 100,658 circuit) | plain
       run_all mode:  python3 check_ref_link.py all [BUILD_DIR]   (the three circuits of the paper, 329 lemmas)

For each lemma l of the variant (sharef.lemma_list()):
  1. the SAT reference B is rebuilt exactly as lemmas.py builds it (lemmas.build_lemma on the candidate's cut
     map); the call of sharef.lemma_apply inside build_lemma is intercepted, recording the node count of B
     before the call, the input words passed to it, and the output words it returns;
  2. the SAW fragment F is rebuilt exactly as gen.py builds it (gen.build_fragment: sharef.lemma_apply on
     fresh input words, IV-valued words as constants);
  3. the substitution sigma maps fragment input bit j of free input word n to the literal of bit j of word n
     in B (one literal per bit, so sigma is a function); constants map to themselves;
  4. check: B's nodes after the recorded count are exactly F's nodes, in order, with operands mapped by
     sigma (and fragment node i mapped to B node count0 + i); the IV-valued input words are identical
     constants in both; every output literal of B (one per alternative literal of each output slot) is the
     image of the fragment output of that slot; and B has no further nodes;
  5. if BUILD_DIR is given (the build directory of run.sh), the AIGER file that SAW read for l is
     byte-identical to aiger.xag_to_aig_bytes(F).
Prints one line per lemma and "REF LINK OK: N of N lemmas" at the end.  No SAT solver is run.
With REF_LINK_SELFTEST=1 one gate type of every fragment is flipped; every lemma must then FAIL.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FLAGS = {"plain": {}, "tail": {"VERIFY_TAIL": "1"}, "split": {"VERIFY_TAIL": "1", "VERIFY_SPLIT": "1"}}
PAPER = [("circuits/proofbase/sha256_compress_96273_proofbase.xag", "split"),
         ("circuits/proofbase/sha256_compress_96324_proofbase.xag", "split"),
         ("circuits/sha256_compress_100658.xag", "tail")]


def run_all(build):
    total = ok = 0
    for circ, var in PAPER:
        env = {k: v for k, v in os.environ.items() if k not in ("VERIFY_TAIL", "VERIFY_SPLIT")}
        env.update(FLAGS[var])
        args = [sys.executable, os.path.abspath(__file__), os.path.join(ROOT, circ), var] + ([build] if build else [])
        out = subprocess.run(args, env=env, capture_output=True, text=True)
        sys.stdout.write(out.stdout); sys.stderr.write(out.stderr)
        last = [ln for ln in out.stdout.splitlines() if ln.startswith("REF LINK")]
        if out.returncode == 0 and last and last[-1].startswith("REF LINK OK"):
            n = int(last[-1].split()[3]); total += n; ok += n
        else:
            total += 1
    print("ALL REF LINKS %s: %d lemmas" % ("OK" if ok == total and total else "FAIL", ok))
    return ok == total and total > 0


def check_one(circ, var, build):
    for k in ("VERIFY_TAIL", "VERIFY_SPLIT"):
        if os.environ.get(k, "0") != FLAGS[var].get(k, "0"):
            raise SystemExit("environment flags do not match variant %s (use the 'all' mode)" % var)
    sys.path[:0] = [os.path.join(ROOT, "verifier"), os.path.join(ROOT, "compositional_proof"), HERE]
    import pickle
    import shutil
    import tempfile
    import lemmas
    import sharef
    import gen
    import aiger
    from xag import read_xag
    from cutmap import build_cutmap

    tmp = tempfile.mkdtemp()
    cp = os.path.join(tmp, os.path.basename(circ))
    shutil.copy(circ, cp)                       # lemmas.py caches the cut map next to the circuit
    x = read_xag(open(cp).read())
    cm = build_cutmap(x, None)
    nodeslots = lemmas.prepare(x, cm)

    cap = []
    real_apply = lemmas.lemma_apply

    def spy(ref, l, inw, style):
        n0 = len(ref.nodes)
        res = real_apply(ref, l, inw, style)
        cap.append((ref, l, {k: list(v) for k, v in inw.items()}, n0, [list(w) for w in res], style))
        return res
    lemmas.lemma_apply = spy

    good = 0
    L = sharef.lemma_list()
    for l in L:
        cap.clear()
        cone, ref, info = lemmas.build_lemma(x, cm, nodeslots, l)
        assert len(cap) == 1, "lemma_apply called %d times" % len(cap)
        B, l2, inwB, n0, resB, style = cap[0]
        assert B is ref and l2 == l and style == "A"
        F, ins, free, outs = gen.build_fragment(sharef, l)
        if os.environ.get("REF_LINK_SELFTEST") and F.nodes:      # negative control: flip one gate type
            op, a, b = F.nodes[len(F.nodes) // 2]
            F.nodes[len(F.nodes) // 2] = ("XOR" if op == "AND" else "AND", a, b)
        errs = []
        # sigma on fragment literals
        sig = {0: 0}
        for k, n in enumerate(free):
            for j in range(32):
                sig[2 * (32 * k + 31 - j + 1)] = inwB[n][j]
        for n in ins:
            if sharef.iv_word(n) is not None and inwB[n] != sharef.cword(sharef.iv_word(n)):
                errs.append("IV word %s differs" % (n,))
        baseF = F.n_inputs + 1
        baseB = B.n_inputs + 1

        def m(lit):
            v = lit >> 1
            if v == 0:
                return lit
            if v < baseF:
                return sig[2 * v] ^ (lit & 1)
            return 2 * (baseB + n0 + (v - baseF)) ^ (lit & 1)
        tail = B.nodes[n0:]
        if len(tail) != len(F.nodes):
            errs.append("node count %d vs fragment %d" % (len(tail), len(F.nodes)))
        else:
            for i, (op, a, b) in enumerate(F.nodes):
                if tail[i] != (op, m(a), m(b)):
                    errs.append("node %d differs" % i); break
        # outputs: fragment output word w of slot outs[w]; B outputs = all alternatives of the output slots
        fo = {}
        for w, name in enumerate(outs):
            for j in range(32):
                fo[(name, j)] = m(F.outputs[32 * w + 31 - j])
        if [fo[(tuple(nm), j)] for (nm, j) in info["out_slots"]] != list(B.outputs):
            errs.append("outputs differ")
        if [[m(lit) for lit in fw] for fw in [[F.outputs[32 * w + 31 - j] for j in range(32)] for w in range(len(outs))]] != resB:
            errs.append("returned words differ")
        if build:
            p = os.path.join(build, var, "frag", gen.lname(l) + ".aig")
            if not os.path.exists(p):
                errs.append("AIGER file missing: %s" % p)
            elif open(p, "rb").read() != aiger.xag_to_aig_bytes(F):
                errs.append("AIGER file differs")
        status = "OK" if not errs else "FAIL " + "; ".join(errs)
        good += not errs
        print("%-10s %-6s nodes %6d  free words %d  sigma bits %4d  out literals %4d  %s"
              % (gen.lname(l), var, len(F.nodes), len(free), 32 * len(free), len(B.outputs), status), flush=True)
    shutil.rmtree(tmp, ignore_errors=True)
    tag = os.path.basename(circ)
    print("REF LINK %s: %d of %d lemmas (%s, variant %s)" % ("OK" if good == len(L) else "FAIL", good, len(L), tag, var))
    return good == len(L)


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "all":
        sys.exit(0 if run_all(sys.argv[2] if len(sys.argv) > 2 else None) else 1)
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    sys.exit(0 if check_one(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None) else 1)
