"""Bit-parallel XAG simulator (Python big ints as lns_) and SPEC mapping tools.

A *ln_* is a Python int whose bit k is the value of a signal on test vector k.
simulate(x, lns_, n) evaluates all n vectors at once.

Interface maps (SPEC <-> circuit) are explicit:
  in_map[c]  = SPEC input index feeding circuit input c, or ("const", v)
  out_map[s] = (circuit output index, complement) giving SPEC output bit s

RANDOM TESTS ARE SMOKE TESTS ONLY.  They never establish equivalence; only a
checked proof (equiv.py) does.
"""

import json
import random
import sys
import time

import sha256_ref as ref
from xag import AND, XAG


def simulate(x: XAG, in_lanes, n):
    if len(in_lanes) != x.n_inputs:
        raise ValueError("wrong number of input lns_")
    mask = (1 << n) - 1
    val = [0] + [v & mask for v in in_lanes]
    for op, a, b in x.nodes:
        va = val[a >> 1]
        if a & 1:
            va ^= mask
        vb = val[b >> 1]
        if b & 1:
            vb ^= mask
        val.append(va & vb if op == AND else va ^ vb)
    return [val[o >> 1] ^ (mask if o & 1 else 0) for o in x.outputs]


def eval_bits(x: XAG, bits):
    """Single-vector evaluation on a list of 0/1 input bits."""
    return [v & 1 for v in simulate(x, [b & 1 for b in bits], 1)]


# --------------------------------------------------------------------------- #
# transposition between lns_ and per-vector integers
# --------------------------------------------------------------------------- #
def lanes_to_vectors(lns_, n):
    """lns_[j] bit k -> vectors[k] bit j."""
    if not lns_:
        return [0] * n
    strs = [format(l, f"0{n}b") for l in lns_]          # char n-1-k is bit k
    rows = zip(*strs)                                     # rows[i] -> vector n-1-i
    vecs = [int("".join(r)[::-1], 2) for r in rows]
    vecs.reverse()
    return vecs


def vectors_to_lanes(vecs, width):
    n = len(vecs)
    strs = [format(v, f"0{width}b") for v in vecs]       # char width-1-j is bit j
    cols = zip(*strs)                                     # cols[i] -> bit width-1-i
    lns_ = [int("".join(c)[::-1], 2) for c in cols]
    lns_.reverse()
    return lns_


# --------------------------------------------------------------------------- #
# SPEC reference functions on packed integers (bit i = SPEC bit i)
# --------------------------------------------------------------------------- #
def _rev(v, nbits):
    return int(format(v, f"0{nbits}b")[::-1], 2)


def _int_to_words(v, nwords):
    """Packed SPEC int (bit i = SPEC bit i) -> u32 words.  SPEC bit i is the
    (i)-th most significant bit of the big-endian concatenation of the words,
    so the packed int is that concatenation bit-reversed."""
    be = _rev(v, 32 * nwords)
    return [(be >> (32 * (nwords - 1 - k))) & ref.MASK for k in range(nwords)]


def _words_to_int(words):
    be = 0
    for w in words:
        be = (be << 32) | w
    return _rev(be, 32 * len(words))


def spec_primary_int(v):
    return _words_to_int(ref.compress_fixed_iv(_int_to_words(v, 16)))


def spec_secondary_int(v):
    w = _int_to_words(v, 24)
    return _words_to_int(ref.compress(w[:8], w[8:]))


def spec_sha256_msg256_int(v):
    """SHA-256 of a 256-bit (32-byte) message, i.e. compress(IV, pad(m))."""
    w = _int_to_words(v, 8)
    blk = ref.pad(b"".join(x.to_bytes(4, "big") for x in w))
    assert len(blk) == 1
    return _words_to_int(ref.compress_fixed_iv(blk[0]))


def spec_primary_noff_int(v):
    """Fixed-IV compression WITHOUT feed-forward (diagnostic hypothesis)."""
    w = _int_to_words(v, 16)
    out = ref.compress_fixed_iv(w)
    return _words_to_int([(o - h) & ref.MASK for o, h in zip(out, ref.IV)])


FUNCTIONS = {
    # name: (n_spec_inputs, function on packed int)
    "primary": (512, spec_primary_int),
    "secondary": (768, spec_secondary_int),
    "sha256_msg256": (256, spec_sha256_msg256_int),
    "primary_no_feedforward": (512, spec_primary_noff_int),
}


# --------------------------------------------------------------------------- #
# candidate interface maps
# --------------------------------------------------------------------------- #
# A word of 32 SPEC bits (SPEC bit j = bit 31-j of the word, MSB first) is laid
# out in the circuit by one of these within-word orders: pos[j] gives the
# circuit offset of SPEC bit j.
WORD_BIT_ORDERS = {
    "msb_first": [j for j in range(32)],
    "lsb_first": [31 - j for j in range(32)],
    # little-endian bytes: byte k of the word (k=0 most significant) placed at 3-k
    "le_bytes_msb_first": [8 * (3 - j // 8) + j % 8 for j in range(32)],
    "le_bytes_lsb_first": [8 * (j // 8) + 7 - j % 8 for j in range(32)],
}


def group_map(n_words, word_order, bit_order):
    """Return pos[s] = circuit offset of SPEC bit s inside a group of n_words."""
    bo = WORD_BIT_ORDERS[bit_order]
    pos = []
    for s in range(32 * n_words):
        w, j = divmod(s, 32)
        cw = w if word_order == "fwd" else n_words - 1 - w
        pos.append(32 * cw + bo[j])
    return pos


def candidate_input_maps(n_circ, n_spec):
    """Yield (description, in_map) for plausible layouts.

    The SPEC input is split into groups (secondary: CV 8 words, block 16 words;
    otherwise one group).  Circuit inputs are the groups in some order, each
    laid out by (word order, within-word bit order); the same layout is used for
    every group.
    """
    if n_circ != n_spec:
        return
    groups = [("cv", 0, 8), ("block", 256, 16)] if n_spec == 768 else [("all", 0, n_spec // 32)]
    orders = [groups, groups[::-1]] if len(groups) > 1 else [groups]
    for gorder in orders:
        for wo in ("fwd", "rev"):
            for bo in WORD_BIT_ORDERS:
                in_map = [None] * n_circ
                off = 0
                for name, sbase, nw in gorder:
                    for s, p in enumerate(group_map(nw, wo, bo)):
                        in_map[off + p] = sbase + s
                    off += 32 * nw
                desc = {"group_order": [g[0] for g in gorder], "word_order": wo,
                        "bit_order": bo}
                yield desc, in_map


def describe_output_map(out_map, n_words=8):
    """Name a discovered output permutation if it matches a standard layout."""
    if any(c for _, c in out_map):
        return None
    perm = [p for p, _ in out_map]
    for wo in ("fwd", "rev"):
        for bo in WORD_BIT_ORDERS:
            if group_map(n_words, wo, bo) == perm:
                return {"word_order": wo, "bit_order": bo}
    return None


# --------------------------------------------------------------------------- #
# evaluation under a map
# --------------------------------------------------------------------------- #
def random_spec_inputs(n_spec, n, seed):
    rng = random.Random(seed)
    return [rng.getrandbits(n_spec) for _ in range(n)]


def eval_under_map(x, in_map, out_map, spec_vecs, n_spec):
    """Evaluate circuit on SPEC-ordered vectors; return SPEC-ordered outputs."""
    n = len(spec_vecs)
    spec_lanes = vectors_to_lanes(spec_vecs, n_spec)
    mask = (1 << n) - 1
    circ_lanes = []
    for m in in_map:
        if isinstance(m, int):
            circ_lanes.append(spec_lanes[m])
        else:
            circ_lanes.append(mask if m[1] else 0)
    outs = simulate(x, circ_lanes, n)
    spec_out = [outs[p] ^ (mask if c else 0) for p, c in out_map]
    return lanes_to_vectors(spec_out, n)


def discover_output_map(circ_out_lanes, ref_out_lanes, mask):
    """Match each SPEC output ln_ to a unique circuit output ln_ (or its
    complement).  Returns out_map or None."""
    idx = {}
    for c, l in enumerate(circ_out_lanes):
        idx.setdefault(l, []).append((c, 0))
        idx.setdefault(l ^ mask, []).append((c, 1))
    out_map, used = [], set()
    for l in ref_out_lanes:
        hits = [h for h in idx.get(l, []) if h[0] not in used]
        if len(hits) != 1:
            return None
        used.add(hits[0][0])
        out_map.append(hits[0])
    return out_map


def discover(x: XAG, functions=None, n=64, seed=1):
    """Try every (function hypothesis, candidate input map); for each, derive the
    output map by ln_ matching.  Returns list of hits."""
    hits = []
    for fname, (n_spec, f) in FUNCTIONS.items():
        if functions and fname not in functions:
            continue
        vecs = random_spec_inputs(n_spec, n, seed)
        ref_lanes = vectors_to_lanes([f(v) for v in vecs], 256)
        spec_lanes = vectors_to_lanes(vecs, n_spec)
        mask = (1 << n) - 1
        for desc, in_map in candidate_input_maps(x.n_inputs, n_spec):
            outs = simulate(x, [spec_lanes[m] for m in in_map], n)
            om = discover_output_map(outs, ref_lanes, mask)
            if om is not None:
                hits.append({"function": fname, "input_layout": desc,
                             "in_map": in_map, "out_map": om,
                             "output_layout": describe_output_map(om)})
    return hits


def known_vectors(fname):
    """Known-answer inputs (packed SPEC ints) for a function."""
    if fname in ("primary", "primary_no_feedforward"):
        # padded "abc" and padded "" (FIPS 180-4 / NIST examples), all-zero, all-one
        blocks = [ref.pad(b"abc")[0], ref.pad(b"")[0], [0] * 16, [ref.MASK] * 16]
        return [_words_to_int(b) for b in blocks]
    if fname == "secondary":
        b1, b2 = ref.pad(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")
        h1 = ref.compress(ref.IV, b1)
        return [_words_to_int(list(ref.IV) + b1), _words_to_int(h1 + b2),
                _words_to_int([0] * 24), _words_to_int([ref.MASK] * 24)]
    if fname == "sha256_msg256":
        return [_words_to_int([0] * 8), _words_to_int([ref.MASK] * 8),
                _words_to_int([int.from_bytes(b"abcd" * 1, "big")] * 8)]
    return []


def verify(x, fname, in_map, out_map, n=100_000, seed=12345, chunk=10_000):
    """Smoke test: compare circuit (under map) with the reference on n random
    vectors plus known vectors.  Returns a dict; mismatches -> 'fail'."""
    n_spec, f = FUNCTIONS[fname]
    t0 = time.time()
    tested, bad = 0, []
    vec_sets = [known_vectors(fname)]
    rng = random.Random(seed)
    left = n
    while left > 0:
        k = min(chunk, left)
        vec_sets.append([rng.getrandbits(n_spec) for _ in range(k)])
        left -= k
    for vecs in vec_sets:
        got = eval_under_map(x, in_map, out_map, vecs, n_spec)
        for v, g in zip(vecs, got):
            if g != f(v):
                bad.append(hex(v))
        tested += len(vecs)
    return {"function": fname, "vectors_tested": tested, "random_vectors": n,
            "known_vectors": len(vec_sets[0]), "seed": seed,
            "mismatches": len(bad), "first_mismatch": bad[:1],
            "status": "pass" if not bad else "fail",
            "note": "random/known-vector simulation is a smoke test, not a proof",
            "seconds": round(time.time() - t0, 1)}


def normalize(x: XAG, in_map, out_map, n_spec) -> XAG:
    """Rebuild x as a circuit in SPEC interface order: input i = SPEC input i,
    output s = SPEC output s.  Unmapped circuit inputs must be constants."""
    y = XAG(n_spec)
    lit = [0]
    for m in in_map:
        lit.append(y.input(m) if isinstance(m, int) else int(m[1]))
    base = x.n_inputs + 1
    for op, a, b in x.nodes:
        la = lit[a >> 1] ^ (a & 1)
        lb = lit[b >> 1] ^ (b & 1)
        lit.append(y._add(op, la, lb))
    y.outputs = []
    for p, c in out_map:
        o = x.outputs[p]
        y.outputs.append(lit[o >> 1] ^ (o & 1) ^ c)
    assert len(lit) == base + len(x.nodes)
    y.validate()
    return y


if __name__ == "__main__":
    from xag import load
    rc = load(sys.argv[1])
    x = rc.to_xag()
    for h in discover(x):
        print(json.dumps({k: v for k, v in h.items() if k not in ("in_map", "out_map")}))
