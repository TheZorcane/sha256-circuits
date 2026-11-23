#!/usr/bin/env python3
"""Soundness: every node of NEW gets a formal affine form over GF(2) in the basis
{1, inputs, A_k (REF AND gates)}. XOR: sum. AND(P,Q) with P,Q in span(F_j,G_j)+const and
a1b2+a2b1=1: PQ = FG + (a1a2+a1c2+a2c1)F + (b1b2+b1c2+b2c1)G + c1c2 (uses x*x=x, valid for
Boolean functions), FG := A_j = value of REF AND j. By induction in NEW topological order,
the form with A_j := function(REF AND j) equals the function of each NEW node (REF's own
forms F_j,G_j are exact under the same substitution, by the REF pass). Equal output forms
=> equal functions. Degenerate planes (F,G dependent mod const, zero) are never indexed =>
fail conservatively. Any matching j is sound (identity holds for any j). Exact big-int
comparison; no hashing."""
import sys, time
from xag2 import parse


def main(refp, newp):
    t0 = time.time()
    n, rn, ro = parse(refp); n2, nn, no = parse(newp)
    if n != n2 or len(ro) != len(no):
        print("FAIL: interface"); return 1
    fm = [0] + [1 << v for v in range(1, n + 1)]
    lf = lambda l: fm[l >> 1] ^ (l & 1)
    ands = []
    for k, a, b in rn:
        if k == 'X':
            fm.append(lf(a) ^ lf(b))
        else:
            F, G = lf(a), lf(b); idx = n + 1 + len(ands); ands.append((F, G)); fm.append(1 << idx)
    rout = [lf(o) for o in ro]
    planes = {}
    for j, (F, G) in enumerate(ands):
        f, g = F & ~1, G & ~1
        if f == 0 or g == 0 or f == g:
            continue
        planes.setdefault(frozenset((f, g, f ^ g)), []).append(j)
    byelem = {}
    for key, js in planes.items():
        for e in key:
            byelem.setdefault(e, []).append(js)
    fn = [0] + [1 << v for v in range(1, n + 1)]
    lg = lambda l: fn[l >> 1] ^ (l & 1)
    rechosen = 0; nA = nX = 0; amb = 0
    for k, a, b in nn:
        P, Q = lg(a), lg(b)
        if k == 'X':
            fn.append(P ^ Q); nX += 1; continue
        nA += 1; res = None
        for js in byelem.get(P & ~1, []):
            for j in js:
                F, G = ands[j]; f, g = F & ~1, G & ~1
                if (Q & ~1) not in (f, g, f ^ g):
                    continue

                def co(X):
                    # X = a*F + b*G + c  (full forms incl. constant bits)
                    for a_ in (0, 1):
                        for b_ in (0, 1):
                            r = X ^ (F if a_ else 0) ^ (G if b_ else 0)
                            if r in (0, 1):
                                return a_, b_, r
                    return None
                cp, cq = co(P), co(Q)
                if not cp or not cq:
                    continue
                a1, b1, c1 = cp; a2, b2, c2 = cq
                if (a1 * b2 + a2 * b1) & 1 != 1:
                    continue
                r = (1 << (n + 1 + j)) ^ (F if (a1 * a2 + a1 * c2 + a2 * c1) & 1 else 0) \
                    ^ (G if (b1 * b2 + b1 * c2 + b2 * c1) & 1 else 0) ^ ((c1 & c2) & 1)
                res = (r, j, len(js)); break
            if res:
                break
        if res is None:
            print("FAIL: NEW AND node #%d has no matching REF plane" % (len(fn) - n - 1)); return 1
        fn.append(res[0]); amb += res[2] > 1
        F, G = ands[res[1]]
        if {P & ~1, Q & ~1} != {F & ~1, G & ~1}:
            rechosen += 1
    bad = [i for i, (x, y) in enumerate(zip(no, ro)) if lg(x) != rout[i]]
    print("REF: AND=%d XOR=%d | NEW: AND=%d XOR=%d | rechosen AND=%d | ambiguous-plane matches=%d | %.1fs"
          % (len(ands), len(rn) - len(ands), nA, nX, rechosen, amb, time.time() - t0))
    if bad:
        print("FAIL: outputs differ in form:", bad[:10], len(bad)); return 1
    print("PASS"); return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
