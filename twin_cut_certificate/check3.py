#!/usr/bin/env python3
"""check3.py REF.xag NEW.xag -- third, independently written GF(2)/ANF equivalence checker (twin-cut certificate).
Written from the spec only (not from lincheck.py / check2.py code).

Representation: each signal is a python int bitmask over the basis {const, x_1..x_n, A_1..A_m}; bit 0 = constant 1,
bit i = x_i, bit n+j = value of the j-th AND gate of REF.  Everything XOR / complement is linear.
REF pass: standard.  Each REF AND j has operand forms F,G; we remember it under the key of each of the three nonzero
constant-free members of span{F',G'} -- as a SET of candidate gates indexed by frozenset plane.
NEW pass: an AND gate with operand forms P,Q: search any REF AND j with P',Q' in the plane of j and P',Q' distinct
nonzero.  Then P = a1 F + b1 G + c1, Q = a2 F + b2 G + c2 (coefficients recovered by exact equality, constants
included).  The product as a polynomial in the *bits* f,g (f,g in {0,1}): evaluate P*Q on the four points
(f,g) and take the ANF by Moebius transform (no closed-form identity used): P*Q = k_fg f g + k_f f + k_g g + k_1.
f*g := A_j (the REF AND gate output; for REF's own operands this is exact).  If k_fg=1 the node form is
k_fg A_j + k_f F + k_g G + k_1; if k_fg=0 the product is affine in (F,G) and the form is k_f F + k_g G + k_1
(sound as well).  Soundness: by induction in NEW's order, form(v) read with A_j := f_j*g_j equals NEW's node.
Outputs compared as full integers (including the constant bit, i.e. complements).
"""
import sys


def read(path):
    toks = open(path).read().split()
    assert toks[:2] == ["xag", "1"] and toks[2] == "inputs" and toks[4] == "nodes"
    n, m = int(toks[3]), int(toks[5])
    p = 6
    nodes = []
    for k in range(m):
        op, a, b = toks[p], int(toks[p + 1]), int(toks[p + 2]); p += 3
        assert op in ("A", "X")
        assert 0 <= a < 2 * (n + 1 + k) and 0 <= b < 2 * (n + 1 + k)
        nodes.append((op == "A", a, b))
    assert toks[p] == "outputs"
    k = int(toks[p + 1]); outs = list(map(int, toks[p + 2:]))
    assert len(outs) == k
    assert all(0 <= o < 2 * (n + 1 + m) for o in outs)
    return n, nodes, outs


def main(rp, np_):
    n, rn, ro = read(rp); n2, nn, no = read(np_)
    assert n == n2 and len(ro) == len(no)
    sig = [0] + [1 << (i + 1) for i in range(n)]
    L = lambda s, l: s[l >> 1] ^ (l & 1)
    ref_ands = []
    for is_and, a, b in rn:
        if is_and:
            F, G = L(sig, a), L(sig, b)
            sig.append(1 << (n + 1 + len(ref_ands)))
            ref_ands.append((F, G))
        else:
            sig.append(L(sig, a) ^ L(sig, b))
    rout = [L(sig, o) for o in ro]
    planes = {}
    for j, (F, G) in enumerate(ref_ands):
        f, g = F & ~1, G & ~1
        if f and g and f != g:
            planes.setdefault(frozenset((f, g, f ^ g)), []).append(j)
    s = [0] + [1 << (i + 1) for i in range(n)]
    nA = 0
    for is_and, a, b in nn:
        P, Q = L(s, a), L(s, b)
        if not is_and:
            s.append(P ^ Q); continue
        pl = frozenset((P & ~1, Q & ~1, (P ^ Q) & ~1))
        cands = planes.get(pl) if len(pl) == 3 and 0 not in pl else None
        if not cands:
            print("FAIL: AND without matching REF plane"); return 1
        j = cands[0]; F, G = ref_ands[j]
        # coefficient recovery (exact)
        def coef(X):
            for a_ in (0, 1):
                for b_ in (0, 1):
                    r = X ^ (F if a_ else 0) ^ (G if b_ else 0)
                    if r in (0, 1) and (a_ or b_):
                        return a_, b_, r
            return None
        cp, cq = coef(P), coef(Q)
        if cp is None or cq is None:
            print("FAIL: coefficient recovery"); return 1
        # truth table of P*Q over (f,g)
        tt = {}
        for f in (0, 1):
            for g in (0, 1):
                pv = (cp[0] & f) ^ (cp[1] & g) ^ cp[2]
                qv = (cq[0] & f) ^ (cq[1] & g) ^ cq[2]
                tt[(f, g)] = pv & qv
        k1 = tt[(0, 0)]
        kf = tt[(1, 0)] ^ k1
        kg = tt[(0, 1)] ^ k1
        kfg = tt[(1, 1)] ^ kf ^ kg ^ k1
        v = (1 << (n + 1 + j)) if kfg else 0
        if kf: v ^= F
        if kg: v ^= G
        v ^= k1
        s.append(v); nA += 1
    nout = [L(s, o) for o in no]
    bad = [i for i in range(len(ro)) if rout[i] != nout[i]]
    if bad:
        print(f"FAIL: {len(bad)} outputs differ, first {bad[0]}"); return 1
    print(f"PASS3: 256-output forms identical; NEW ANDs {nA}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
