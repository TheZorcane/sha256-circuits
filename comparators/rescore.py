"""Helpers for re-scoring circuits (library module, used by textbook.py).

Uses verifier/{xag,cost,sim}.py read-only for parsing, canonical counts and
SHA-256 simulation; own refs.py for SHA-1 / SHA-512.  not_tally() computes a
NOT-counting tally (NIST convention: AND, XOR, XNOR, NOT all count 1) for any
XAG, so XAG circuits can be stated in the published NIST convention.
"""
import hashlib
import json
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
sys.path.insert(0, HERE)
import cost  # noqa: E402
import sim  # noqa: E402
from xag import load, read_xag, AND, XOR  # noqa: E402
import refs  # noqa: E402

SEED = 20261119
NSIM = int(os.environ.get("NSIM", "20000"))


def sha256f(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def not_tally(x):
    """Gates in basis {AND, XOR, XNOR, NOT} (each counts 1) needed to realise
    XAG x without restructuring: XOR nodes may be emitted as XNOR (free choice of
    polarity, operand complements absorbed); a NOT is needed once per signal for
    which some AND operand / primary output needs a polarity the signal does not
    natively provide."""
    base = x.n_inputs + 1
    need = {}  # var -> set of polarities required by AND operands / outputs
    for op, a, b in x.nodes:
        if op == AND:
            for l in (a, b):
                if l >> 1:
                    need.setdefault(l >> 1, set()).add(l & 1)
    for o in x.outputs:
        if o >> 1:
            need.setdefault(o >> 1, set()).add(o & 1)
    nots = 0
    for v, pols in need.items():
        is_xor = v >= base and x.nodes[v - base][0] == XOR
        if is_xor:
            nots += 1 if len(pols) == 2 else 0
        else:
            nots += 1 if 1 in pols else 0
    c = x.counts()
    return {"AND": c["AND"], "XOR_XNOR": c["XOR"], "NOT": nots,
            "tally_total": c["total"] + nots}


# ---------------- generic interface discovery for SHA-1 / SHA-512 ----------
def word_bits(words, w):
    v = 0
    for i, x in enumerate(words):
        for j in range(w):
            if (x >> (w - 1 - j)) & 1:
                v |= 1 << (i * w + j)
    return v


def bits_word(v, nw, w):
    out = []
    for i in range(nw):
        x = 0
        for j in range(w):
            x = (x << 1) | ((v >> (i * w + j)) & 1)
        out.append(x)
    return out


def make_fn(kind):
    if kind == "sha1_fix":
        return 512, 32, lambda v: word_bits(refs.sha1_compress(refs.SHA1_IV, bits_word(v, 16, 32)), 32), 160, [("all", 0, 16)]
    if kind == "sha512_var":
        def f(v):
            ws = bits_word(v, 24, 64)
            return word_bits(refs.sha512_compress(ws[:8], ws[8:]), 64)
        return 1536, 64, f, 512, [("cv", 0, 8), ("block", 512, 16)]
    if kind == "sha512_fix":
        return 1024, 64, lambda v: word_bits(refs.sha512_compress(refs.SHA512_IV, bits_word(v, 16, 64)), 64), 512, [("all", 0, 16)]
    raise KeyError(kind)


def cand_maps(n, w, groups):
    orders = [groups, groups[::-1]] if len(groups) > 1 else [groups]
    for go in orders:
        for wo in ("fwd", "rev"):
            for bo in ("msb", "lsb"):
                in_map = [None] * n
                off = 0
                for name, sb, nw in go:
                    for s in range(nw * w):
                        wd, j = divmod(s, w)
                        cw = wd if wo == "fwd" else nw - 1 - wd
                        cj = j if bo == "msb" else w - 1 - j
                        in_map[off + cw * w + cj] = sb + s
                    off += nw * w
                yield f"{[g[0] for g in go]}/{wo}/{bo}", in_map


def verify_generic(x, kind, n=NSIM, seed=SEED):
    n_spec, w, f, n_out, groups = make_fn(kind)
    rng = random.Random(seed)
    probe = [rng.getrandbits(n_spec) for _ in range(64)]
    ref_l = sim.vectors_to_lanes([f(v) for v in probe], n_out)
    spec_l = sim.vectors_to_lanes(probe, n_spec)
    mask = (1 << 64) - 1
    hits = []
    for desc, im in cand_maps(x.n_inputs, w, groups):
        outs = sim.simulate(x, [spec_l[m] for m in im], 64)
        om = sim.discover_output_map(outs, ref_l, mask)
        if om is not None:
            hits.append((desc, im, om))
    if len(hits) != 1:
        return {"status": "no unique interface", "hits": len(hits)}
    desc, im, om = hits[0]
    bad, done, t0 = 0, 0, time.time()
    while done < n:
        k = min(2000, n - done)
        vecs = [rng.getrandbits(n_spec) for _ in range(k)]
        if done == 0:
            vecs[:3] = [0, (1 << n_spec) - 1, int("5a" * (n_spec // 8), 16)]
        lns_ = sim.vectors_to_lanes(vecs, n_spec)
        outs = sim.simulate(x, [lns_[m] for m in im], k)
        km = (1 << k) - 1
        got = sim.lanes_to_vectors([outs[p] ^ (km if c else 0) for p, c in om], k)
        bad += sum(1 for v, g in zip(vecs, got) if g != f(v))
        done += k
    return {"status": "pass" if bad == 0 else "FAIL", "interface": desc,
            "vectors": done, "mismatches": bad, "seed": seed,
            "out_map_identity": all(p == i and c == 0 for i, (p, c) in enumerate(om)),
            "seconds": round(time.time() - t0, 1)}


def verify_sha256(x, fname, n=NSIM):
    hits = sim.discover(x, functions=[fname])
    if len(hits) != 1:
        return {"status": "no unique interface", "hits": len(hits)}
    h = hits[0]
    r = sim.verify(x, fname, h["in_map"], h["out_map"], n=n, seed=SEED)
    r["interface"] = h["input_layout"]
    return r
