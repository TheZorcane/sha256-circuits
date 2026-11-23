"""Canonicalization and cost report for XAGs.

canonicalize(x) applies, in one topological pass:
  * NOT absorption (already implicit: complements live on edges);
  * XOR complement normalisation: XOR(~a, b) = ~XOR(a, b), so XOR operands are
    stored uncomplemented and the complement moves to the output edge;
  * constant propagation and trivial identities
      AND(a,0)=0  AND(a,1)=a  AND(a,a)=a  AND(a,~a)=0
      XOR(a,0)=a  XOR(a,1)=~a XOR(a,a)=0  XOR(a,~a)=1;
  * operand ordering (a <= b) and structural hashing (identical (op,a,b)
    nodes merged);
then removes dead nodes and renumbers the survivors in a deterministic DFS
post-order from the outputs (output order; operands visited in order of a
bottom-up structural digest, not source order), so that the serialization --
and its sha256 -- depends neither on the gate order nor on the operand order
in the source file.  (Two circuits that are identical up to wire renaming,
gate reordering and operand swaps get the same hash; the converse holds
trivially since the serialization determines the circuit.)  The pass is
idempotent (tested).

No other rewriting is done: the canonical count is *not* an optimised count.
"""

import hashlib
import sys

from xag import AND, XOR, XAG, RawCircuit, write_xag


def _simplify(x: XAG) -> XAG:
    y = XAG(x.n_inputs)
    m = list(range(0, 2 * (x.n_inputs + 1), 2))   # old var -> new lit
    strash = {}

    def tr(l):
        return m[l >> 1] ^ (l & 1)

    for op, a, b in x.nodes:
        a, b = tr(a), tr(b)
        neg = 0
        if op == XOR:
            neg = (a & 1) ^ (b & 1)
            a, b = a & ~1, b & ~1
            if a > b:
                a, b = b, a
            if a == 0:
                r = b
            elif a == b:
                r = 0
            else:
                r = None
        else:
            if a > b:
                a, b = b, a
            if a == 0:
                r = 0
            elif a == 1:
                r = b
            elif a == b:
                r = a
            elif a == b ^ 1:
                r = 0
            else:
                r = None
        if r is None:
            key = (op, a, b)
            r = strash.get(key)
            if r is None:
                r = y._add(op, a, b)
                strash[key] = r
        m.append(r ^ neg)
    y.outputs = [tr(o) for o in x.outputs]
    return y


def _merkle(x: XAG):
    """Per-var structural digest, invariant under operand order and gate order:
    d(input i) = H('i', i); d(node) = H(op, sorted(d(a)+neg_a, d(b)+neg_b))."""
    def H(*parts):
        return hashlib.blake2b(repr(parts).encode(), digest_size=16).digest()
    d = [H("c")] + [H("i", i) for i in range(x.n_inputs)]
    for op, a, b in x.nodes:
        ka, kb = d[a >> 1] + bytes([a & 1]), d[b >> 1] + bytes([b & 1])
        d.append(H(op, *sorted((ka, kb))))
    return d


def _dfs_renumber(x: XAG) -> XAG:
    base = x.n_inputs + 1
    dig = _merkle(x)
    new_var = {v: v for v in range(base)}
    y = XAG(x.n_inputs)
    for o in x.outputs:
        stack = [(o >> 1, False)]
        while stack:
            v, expanded = stack.pop()
            if v in new_var:
                continue
            op, a, b = x.nodes[v - base]
            if expanded:
                na = 2 * new_var[a >> 1] | (a & 1)
                nb = 2 * new_var[b >> 1] | (b & 1)
                if na > nb:
                    na, nb = nb, na
                new_var[v] = y._add(op, na, nb) >> 1
            else:
                # visit operands in digest order (not source order) so the
                # numbering is independent of how the file wrote the operands
                ka, kb = dig[a >> 1] + bytes([a & 1]), dig[b >> 1] + bytes([b & 1])
                first, second = (a, b) if ka <= kb else (b, a)
                stack.append((v, True))
                stack.append((second >> 1, False))
                stack.append((first >> 1, False))
    y.outputs = [2 * new_var[o >> 1] | (o & 1) for o in x.outputs]
    return y


def canonicalize(x: XAG) -> XAG:
    x.validate()
    y = _dfs_renumber(_simplify(x))
    y.validate()
    return y


def serialize(x: XAG) -> bytes:
    return write_xag(x).encode()


def xag_hash(x: XAG) -> str:
    return hashlib.sha256(serialize(x)).hexdigest()


def depth_stats(x: XAG):
    base = x.n_inputs + 1
    depth = [0] * x.n_vars()
    adepth = [0] * x.n_vars()
    fan = [0] * x.n_vars()
    for k, (op, a, b) in enumerate(x.nodes):
        va, vb = a >> 1, b >> 1
        fan[va] += 1
        fan[vb] += 1
        depth[base + k] = 1 + max(depth[va], depth[vb])
        adepth[base + k] = (op == AND) + max(adepth[va], adepth[vb])
    for o in x.outputs:
        fan[o >> 1] += 1
    outs = [o >> 1 for o in x.outputs]
    return {
        "depth": max((depth[v] for v in outs), default=0),
        "and_depth": max((adepth[v] for v in outs), default=0),
        # fanout = number of node operand slots + output slots referencing a var
        "max_fanout_node": max(fan[base:], default=0),
        "max_fanout_input": max(fan[1:base], default=0),
    }


def report_xag(x: XAG) -> dict:
    r = dict(x.counts())
    r.update(depth_stats(x))
    s = serialize(x)
    r["serialized_bytes"] = len(s)
    r["sha256"] = hashlib.sha256(s).hexdigest()
    return r


def report(rc: RawCircuit) -> dict:
    """Full report: raw (file) counts, direct XAG translation, canonical XAG."""
    direct = rc.to_xag()
    canon = canonicalize(direct)
    return {
        "raw_counts_by_type": rc.raw_counts(),
        "n_inputs": len(rc.inputs),
        "n_outputs": len(rc.outputs),
        "direct_xag": report_xag(direct),
        "canonical_xag": report_xag(canon),
    }


if __name__ == "__main__":
    import json
    from xag import load, read_xag
    p = sys.argv[1]
    if p.endswith(".xag"):
        x = read_xag(open(p).read())
        print(json.dumps({"input_xag": report_xag(x),
                          "canonical_xag": report_xag(canonicalize(x))}, indent=2))
    else:
        print(json.dumps(report(load(p)), indent=2))
