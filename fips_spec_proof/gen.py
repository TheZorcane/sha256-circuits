"""Generate the inputs of the SAW proofs that link the verifier's reference to the Cryptol FIPS 180-4 spec.

usage: python3 gen.py BUILD_DIR

For each verifier variant (plain, tail, split) this imports the verifier's own reference generator
../compositional_proof/sharef.py (read-only, with the variant's environment flags) and writes:
  BUILD/<var>/frag/<kind>_<i>.aig    reference fragment circuit = sharef.lemma_apply(..., "A") on fresh input
                                     words (IV-valued input words become constants, as in lemmas.py).  The
                                     SAT reference of lemmas.py consists of the same lemma_apply gates applied
                                     to the cone's input literals (plus, for expanded input words, the
                                     candidate's own gates below them); check_ref_link.py checks, for every
                                     lemma, that these gates are node for node this fragment under a
                                     consistent input substitution and that the .aig file SAW reads is
                                     byte-identical to the converted fragment
  BUILD/<var>/Lemmas_<var>_S.cry     lem_<kind>_<i>: Cryptol function of each lemma (Frag.cry), and
                                     composed : [16][32] -> [8][32], the composition of the lem_* functions
                                     over sharef's wiring (lemma_list / lemma_in_out / iv_word), generated
  BUILD/<var>/Lemmas_<var>_I.cry     the same with SHA2Internal's word operations (FragI.cry)
  BUILD/<var>/frags_{S,I}.saw        one SAW theorem per lemma: aig == lem_*
plus BUILD/FragI.cry, the negative-control mutant and BUILD/gen_info.json.  Every AIGER file is validated
by random simulation against verifier/sim.py.

Variants (environment flags of sharef.py):
  plain  {}                                   120 lemmas: schedule, rounds 0..63, 8 feed-forward words
  tail   VERIFY_TAIL=1                        109 lemmas: rounds 62-63 + feed-forward fused (one tail lemma)
  split  VERIFY_TAIL=1 VERIFY_SPLIT=1         110 lemmas: the fused tail split in two lemmas

Bit/word mapping (verifier/SPEC.md <-> Cryptol): a sharef word is a list of 32 literals, LSB first.
Fragment AIGER input 32*k + (31-j) = bit j of the k-th free input word; outputs likewise.  SAW's read_aig
maps AIGER input i to element i of a Cryptol [n] (element 0 = most significant), so the AIGER term takes
`join [w_0, w_1, ...]` with w_k : [32] an ordinary big-endian Cryptol word.  This is SPEC.md's order.
"""
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SHAREF = os.path.join(ROOT, "compositional_proof", "sharef.py")
sys.path.insert(0, os.path.join(ROOT, "verifier"))
sys.path.insert(0, HERE)
import aiger                      # noqa: E402
import cost                       # noqa: E402  (hashing only)
import sha256_ref as S            # noqa: E402
from xag import XAG               # noqa: E402

FLAGS = ("VERIFY_TAIL", "VERIFY_SPLIT")
VARIANTS = {
    "plain": {},
    "tail": {"VERIFY_TAIL": "1"},
    "split": {"VERIFY_TAIL": "1", "VERIFY_SPLIT": "1"},
}


def load_sharef(env, tag):
    saved = {k: os.environ.pop(k) for k in FLAGS if k in os.environ}
    os.environ.update(env)
    spec = importlib.util.spec_from_file_location(f"sharef_{tag}", SHAREF)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    for k in FLAGS:
        os.environ.pop(k, None)
    os.environ.update(saved)
    return m


def vname(n):
    kind, t = n
    return f"{kind}_{t}" if t >= 0 else f"{kind}_m{-t}"


def lname(l):
    return f"{l[0]}_{l[1]}"


# Cryptol expression of each lemma kind, mirroring sharef.lemma_apply's dispatch.  args: Cryptol expressions
# of the declared input words (lemma_in_out order); result: list of output words (lemma_in_out order).
def cry_apply(kind, i, a):
    if kind == "sched":                       # sched_step(w2, w7, w15, w16)
        return f"[schedStep {a[0]} {a[1]} {a[2]} {a[3]}]"
    if kind == "round":                       # round_step(state8, W, K[i]) -> (a', e')
        return f"(\\(p, q) -> [p, q]) (roundStep [{', '.join(a[:8])}] (K @ ({i} : [8])) {a[8]})"
    if kind == "final":                       # final_add(IV[i], w)
        return f"[finalAdd (IV @ ({i} : [4])) {a[0]}]"
    if kind == "tail":                        # a61..58, e61..58, W60,W55,W47,W46,W61,W56,W48
        return f"tailFn [{', '.join(a[:8])}] {' '.join(a[8:15])}"
    if kind == "r62f":
        return f"(\\(p, q) -> [p, q]) (r62fFn [{', '.join(a[:8])}] {' '.join(a[8:12])})"
    if kind == "tail2":
        return f"tail2Fn [{', '.join(a[:8])}] {' '.join(a[8:12])}"
    raise ValueError(kind)


def iv_expr(sh, n):
    v = sh.iv_word(n)
    if v is None:
        return None
    kind, t = n
    idx = (-t - 1) if kind == "a" else 4 + (-t - 1)
    assert S.IV[idx] == v
    return f"(IV @ ({idx} : [4]))"


def build_fragment(sh, l):
    ins, outs = sh.lemma_in_out(l)
    free = [n for n in ins if sh.iv_word(n) is None]
    x = XAG(32 * len(free))
    inw = {}
    for n in ins:
        v = sh.iv_word(n)
        if v is not None:
            inw[n] = sh.cword(v)                    # IV words are constants, as in lemmas.py
        else:
            k = free.index(n)
            inw[n] = [x.input(32 * k + 31 - j) for j in range(32)]
    res = sh.lemma_apply(x, l, inw, "A")
    assert len(res) == len(outs)
    x.outputs = [w[31 - j] for w in res for j in range(32)]
    return x, ins, free, outs


def concrete_vectors(x, n, seed):
    """(input int, output int) pairs computed by verifier/sim.py on the XAG (not on the AIGER file).
    Checked in SAW against the read_aig term: an end-to-end test of converter, reader and bit order."""
    import random
    import sim
    r = random.Random(seed)
    out = []
    for _ in range(n):
        bits = [r.getrandbits(1) for _ in range(x.n_inputs)]
        ob = sim.eval_bits(x, bits)
        out.append((int("".join(map(str, bits)), 2), int("".join(map(str, ob)), 2)))
    return out


def gen_variant(var, build, style):
    """style "S": Frag.cry (Specification-form word operations); "I": FragI.cry (SHA2Internal's operations)."""
    env = VARIANTS[var]
    sh = load_sharef(env, var)
    d = os.path.join(build, var)
    os.makedirs(os.path.join(d, "frag"), exist_ok=True)
    cry = [f"// GENERATED by gen.py from compositional_proof/sharef.py, flags {env}. Do not edit.",
           f"module Lemmas_{var}_{style} where", "import Frag" + ("I" if style == "I" else ""), ""]
    saw = [f"// GENERATED: fragment theorems for variant {var}, style {style}", f'import "Lemmas_{var}_{style}.cry";', ""]
    comp_lines = []
    man = []
    for l in sh.lemma_list():
        x, ins, free, outs = build_fragment(sh, l)
        b = aiger.xag_to_aig_bytes(x)
        ok = aiger.validate(x, b)
        assert ok, l
        open(os.path.join(d, "frag", lname(l) + ".aig"), "wb").write(b)
        nI, nO = 32 * len(free), 32 * len(outs)
        args = [iv_expr(sh, n) or f"i{free.index(n)}" for n in ins]
        pat = "[" + ", ".join(f"i{k}" for k in range(len(free))) + "]"
        cry += [f"lem_{lname(l)} : [{nI}] -> [{nO}]",
                f"lem_{lname(l)} x = join ({cry_apply(l[0], l[1], args)})",
                f"  where {pat} = split x : [{len(free)}][32]", ""]
        vecs = concrete_vectors(x, 2, seed=len(man))
        rv = f"r_{lname(l)}"
        saw += [f'{rv} <- read_aig "frag/{lname(l)}.aig";',
                f'print "F{style}_{var}_{lname(l)}: aig(sharef.lemma_apply {l}) == lem_{lname(l)}  [{nI} -> {nO} bits]";']
        saw += [f'prove_print rme {{{{ {rv} 0x{vi:0{nI // 4}x} == 0x{vo:0{nO // 4}x} }}}};' for vi, vo in vecs]
        saw += [f'time (prove_print yices {{{{ \\(x : [{nI}]) -> {rv} x == lem_{lname(l)} x }}}});',
                f'print "PROVED F{style}_{var}_{lname(l)}";', ""]
        fargs = "[" + ", ".join(vname(n) for n in free) + "]"
        comp_lines.append(f"    [{', '.join(vname(n) for n in outs)}] = split (lem_{lname(l)} (join {fargs})) : [{len(outs)}][32]")
        man.append({"lemma": lname(l), "inputs": [vname(n) for n in ins], "free": [vname(n) for n in free],
                    "outputs": [vname(n) for n in outs], "nodes": len(x.nodes),
                    "xag_sha256": cost.xag_hash(x), "aig_validated_vs_sim": ok})
    win = ", ".join(f"W_{t}" for t in range(16))
    cry += ["/** Composition of the lemma functions over the verifier's wiring (sharef.lemma_list order). */",
            "composed : [16][32] -> [8][32]",
            f"composed [{win}] = [{', '.join(f'out_{k}' for k in range(8))}]",
            "  where"] + comp_lines + [""]
    open(os.path.join(d, f"Lemmas_{var}_{style}.cry"), "w").write("\n".join(cry))
    saw.append(f'print "ALL {len(man)} FRAGMENT THEOREMS F{style} PROVED FOR VARIANT {var}";')
    open(os.path.join(d, f"frags_{style}.saw"), "w").write("\n".join(saw))
    json.dump(man, open(os.path.join(d, f"fragments_{style}.json"), "w"), indent=0)
    return sh, man


def inline_composition(sh):
    """Rebuild the whole reference by inlining the proved fragment XAGs along the wiring, as
    sharef.full_primary chains lemma_apply; used to check that full_primary("A") is that composition."""
    x = XAG(512)
    vals = {("W", t): [x.input(32 * t + 31 - j) for j in range(32)] for t in range(16)}
    outs = [None] * 8
    for l in sh.lemma_list():
        f, ins, free, onames = build_fragment(sh, l)
        mp = {0: 0}
        for k, n in enumerate(free):
            for j in range(32):
                mp[32 * k + 31 - j + 1] = vals[n][j]

        def m(lit):
            return mp[lit >> 1] ^ (lit & 1)
        for i, (op, a, b) in enumerate(f.nodes):
            mp[f.n_inputs + 1 + i] = x.AND(m(a), m(b)) if op == "AND" else x.XOR(m(a), m(b))
        res = [[m(f.outputs[32 * w + 31 - j]) for j in range(32)] for w in range(len(onames))]
        for n, w in zip(onames, res):
            if n[0] == "out":
                outs[n[1]] = w
            else:
                vals[n] = w
    x.outputs = [outs[k][31 - j] for k in range(8) for j in range(32)]
    return x


def write_aig_checked(x, path):
    b = aiger.xag_to_aig_bytes(x)
    ok = aiger.validate(x, b, nvec=128)
    assert ok, path
    open(path, "wb").write(b)
    return {"file": os.path.basename(path), "nodes": len(x.nodes), "xag_sha256": cost.xag_hash(x),
            "aig_validated_vs_sim": ok}


def make_fragI(build):
    """FragI.cry = Frag.cry with the six word functions replaced by SHA2Internal's exported ones."""
    src = open(os.path.join(HERE, "cry", "Frag.cry")).read()
    head, rest = src.split("ROTR : {n}", 1)
    body = rest[rest.index("/** W_t ="):]
    ops = "\n".join(["Sigma0 = SHAI::S0", "Sigma1 = SHAI::S1", "sigma0 = SHAI::s0", "sigma1 = SHAI::s1",
                     "Ch = SHAI::Ch", "Maj = SHAI::Maj", ""])
    head = head.replace("module Frag where", "module FragI where")
    out = ("// GENERATED by gen.py from Frag.cry: same functions, but the word operations are SHA2Internal's\n"
           "// exported S0/S1/s0/s1/Ch/Maj.  Used for theorem W_I (equality with processBlock_Common).\n"
           + head + ops + "\n" + body)
    open(os.path.join(build, "FragI.cry"), "w").write(out)


def main(build):
    os.makedirs(build, exist_ok=True)
    make_fragI(build)
    info = {}
    shs = {}
    for var in VARIANTS:
        shs[var], man = gen_variant(var, build, "S")
        gen_variant(var, build, "I")
        info[var] = len(man)
    # the verifier's whole reference is gate-for-gate the inlined composition of the proved fragments
    fpA = shs["plain"].full_primary("A")
    info["full_reference"] = {"nodes": len(fpA.nodes), "xag_sha256": cost.xag_hash(fpA)}
    comp = inline_composition(shs["plain"])
    info["full_reference"]["identical_to_inlined_fragments_plain"] = (
        comp.nodes == fpA.nodes and comp.outputs == fpA.outputs)
    # negative control: one gate of the round-10 reference fragment changed from AND to XOR
    x, ins, free, outs = build_fragment(shs["plain"], ("round", 10))
    k = [i for i, n in enumerate(x.nodes) if n[0] == "AND"][len(x.nodes) // 7]
    op, a, b = x.nodes[k]
    x.nodes[k] = ("XOR", a, b)
    info["mutant_round_10"] = write_aig_checked(x, os.path.join(build, "mutant_round_10.aig"))
    info["mutant_round_10"]["mutated_node"] = k
    json.dump(info, open(os.path.join(build, "gen_info.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
