#!/usr/bin/env python3
"""lincheck.py REF.xag NEW.xag  --  independent GF(2) linear-algebra equivalence checker (linear-layer optimiser).

Stand-alone on purpose: own parser, no imports from the optimiser or from verifier/.

Claim checked: NEW computes the same Boolean function as REF on every input.

Method.  In any XAG every literal is an affine form  c + sum_i x_i + sum_j A_j  over GF(2) in the
primary inputs x_i and the outputs A_j of the AND gates (XOR gates and complemented edges are affine).
We represent a form as a python int (bit 0: constant, bit 1+i: input i, bit 1+n+j: REF's AND j).

1. REF: compute the form of every node in REF's own AND basis.  For each REF AND j with operand
   forms F, G (as literals, constants included) record the *plane* key
   sorted{F', G', F'+G'} (primes: constant bit cleared) -> (F, G, basis vector A_j).
2. NEW: process nodes in topological (file) order.  An XOR node gets the XOR of its operand forms.
   An AND node with operand forms P, Q must hit a REF AND j whose plane key equals
   sorted{P', Q', P'+Q'}.  Then P = a1 F + b1 G + c1 and Q = a2 F + b2 G + c2 for unique bits
   (verified by exact vector equality), with a1 b2 + a2 b1 = 1, and since F*F = F, G*G = G for
   Boolean functions
       P*Q = F*G + (a1 a2 + a1 c2 + a2 c1) F + (b1 b2 + b1 c2 + b2 c1) G + c1 c2 .
   The node gets the form  A_j + (..)F + (..)G + c1 c2  (an identity between functions, so the
   invariant below is preserved).  When P, Q are REF's own operands (the usual case) this is A_j.
   If no REF AND has that plane the check FAILS (no claim is made).
3. Every output literal of NEW must have exactly the same form as the same output of REF.
The PASS line also counts the NEW ANDs whose correction (..)F + (..)G + c1 c2 is zero.  For such an AND,
NEW computes exactly the function of the matched REF AND (plane agreement alone does not give this:
F*G could become F*(F+G) = F*G + F).  Output equality does not depend on this count.

Soundness (induction over NEW's file order, which is a topological order - validated here).
Invariant: for every NEW node v, form(v), read as an affine expression in the inputs and in the
*functions computed by REF's AND gates*, equals the function NEW computes at v.  Inputs: trivial.
XOR: affine forms add.  AND with operand forms P, Q: by the invariant the operands compute the
functions denoted by P and Q; REF's AND j computes A_j = F*G where F, G denote the functions of its
operands (step 1 is the invariant for REF itself, a tautology in REF's own basis); the displayed
polynomial identity then gives the function of NEW's AND as an affine form.  Outputs: equal forms
denote equal functions.  The converse does not hold (different forms can denote equal functions), so
the checker is sound but incomplete; it never reports PASS for inequivalent circuits.  No hashing
of forms into shorter digests is used: dictionary keys are the full integers.
The check does not require the AND correspondence to be injective or surjective.
"""
import sys


def parse(path):
    t = open(path).read().split()
    assert t[0] == "xag" and t[1] == "1" and t[2] == "inputs" and t[4] == "nodes"
    n, m = int(t[3]), int(t[5])
    nodes = []
    i = 6
    for k in range(m):
        op, a, b = t[i], int(t[i + 1]), int(t[i + 2]); i += 3
        assert op in ("A", "X")
        lim = 2 * (n + 1 + k)
        assert 0 <= a < lim and 0 <= b < lim, "not topologically ordered"
        nodes.append((op, a, b))
    assert t[i] == "outputs"
    k = int(t[i + 1]); outs = [int(z) for z in t[i + 2:i + 2 + k]]
    assert len(outs) == k and i + 2 + k == len(t)
    assert all(0 <= o < 2 * (n + 1 + m) for o in outs)
    return n, nodes, outs


def main():
    n, rn, ro = parse(sys.argv[1])
    n2, nn, no = parse(sys.argv[2])
    assert n == n2 and len(ro) == len(no)

    def lf(form, l):
        return form[l >> 1] ^ (l & 1)

    def pkey(f, g):
        f0, g0 = f & ~1, g & ~1
        if not f0 or not g0 or f0 == g0:
            return None
        return tuple(sorted((f0, g0, f0 ^ g0)))

    # 1. REF
    form = [0] + [1 << (1 + i) for i in range(n)]
    key = {}
    j = 0
    for op, a, b in rn:
        if op == "X":
            form.append(lf(form, a) ^ lf(form, b))
        else:
            f, g = lf(form, a), lf(form, b)
            k = pkey(f, g)
            if k is not None:
                key.setdefault(k, (f, g, 1 << (1 + n + j)))
            form.append(1 << (1 + n + j)); j += 1
    ref_out = [lf(form, o) for o in ro]
    n_and_ref = j

    def coords(p, f, g):
        for a1, b1 in ((1, 0), (0, 1), (1, 1)):
            r = p ^ (f if a1 else 0) ^ (g if b1 else 0)
            if r in (0, 1):
                return a1, b1, r
        return None

    # 2. NEW
    form = [0] + [1 << (1 + i) for i in range(n)]
    used = set(); n_and_new = 0; replaned = 0; zero_corr = 0
    for idx, (op, a, b) in enumerate(nn):
        if op == "X":
            form.append(lf(form, a) ^ lf(form, b))
            continue
        p, q = lf(form, a), lf(form, b)
        k = pkey(p, q)
        hit = key.get(k) if k is not None else None
        if hit is None:
            print(f"FAIL: AND node {idx} of NEW has an operand plane matching no AND of REF")
            return 1
        f, g, aj = hit
        c1, c2 = coords(p, f, g), coords(q, f, g)
        if c1 is None or c2 is None or (c1[0] & c2[1]) ^ (c2[0] & c1[1]) != 1:
            print(f"FAIL: AND node {idx} of NEW: operands are not an independent pair of the plane")
            return 1
        a1, b1, k1 = c1; a2, b2, k2 = c2
        cf = (a1 & a2) ^ (a1 & k2) ^ (a2 & k1)
        cg = (b1 & b2) ^ (b1 & k2) ^ (b2 & k1)
        corr = (f if cf else 0) ^ (g if cg else 0) ^ (k1 & k2)
        zero_corr += corr == 0
        form.append(aj ^ corr)
        used.add(aj); n_and_new += 1
        if {p & ~1, q & ~1} != {f & ~1, g & ~1}:
            replaned += 1
    # 3. outputs
    new_out = [lf(form, o) for o in no]
    bad = [i for i in range(len(ro)) if ref_out[i] != new_out[i]]
    if bad:
        print(f"FAIL: {len(bad)} outputs differ as forms, first {bad[0]}")
        return 1
    print(f"PASS: outputs {len(ro)} identical as affine forms; REF ANDs {n_and_ref}, NEW ANDs {n_and_new}, "
          f"distinct REF ANDs matched {len(used)}, with re-chosen operand pair {replaned}; "
          f"NEW XORs {len(nn) - n_and_new}; ANDs equal to their matched REF AND (zero correction) "
          f"{zero_corr} of {n_and_new}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
