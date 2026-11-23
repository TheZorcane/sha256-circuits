"""XAG -> binary AIGER (aig 1.9 subset, combinational) and an independent AIGER reader/simulator.

XAG literals use the AIGER convention already (lit = 2*var + neg, var 0 = FALSE, vars 1..n = inputs),
so inputs keep their numbers.  Each XAG node becomes:
    AND(a,b)  -> one AIGER AND gate
    XOR(a,b)  -> 3 AIGER AND gates:  n1 = AND(a,b); n2 = AND(~a,~b); x = AND(~n1,~n2)
                 (x = ~(a&b) & ~(~a&~b) = a XOR b)
Input i of the XAG is AIGER input i (in file order); output k of the XAG is AIGER output k.
This module deliberately does not use verifier/xag.py.

The reader (read_aig / simulate_aig) is a separate, straightforward decoder used only to validate the
writer by simulation against verifier/sim.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VERIFIER = os.path.join(HERE, "..", "verifier")


def _enc(n):
    out = bytearray()
    while n >= 0x80:
        out.append((n & 0x7F) | 0x80)
        n >>= 7
    out.append(n)
    return bytes(out)


def xag_to_aig_bytes(x):
    I = x.n_inputs
    amap = {0: 0}                      # XAG var -> AIGER literal (positive)
    for i in range(1, I + 1):
        amap[i] = 2 * i
    gates = []                         # (lhs, rhs0, rhs1)
    nxt = [I + 1]

    def lit(l):
        return amap[l >> 1] ^ (l & 1)

    strash = {}

    def mk_and(a, b):
        # structural hashing + trivial identities (SAW's reader goes through ABC, which strashes on read;
        # duplicate / trivial gates in the file make SAW fail with "Invalid literal")
        if a < b:
            a, b = b, a
        if b == 0:
            return 0
        if b == 1 or a == b:
            return a
        if a == b ^ 1:
            return 0
        if (a, b) in strash:
            return strash[(a, b)]
        lhs = 2 * nxt[0]
        nxt[0] += 1
        gates.append((lhs, a, b))
        strash[(a, b)] = lhs
        return lhs

    for k, (op, a, b) in enumerate(x.nodes):
        la, lb = lit(a), lit(b)
        if la > lb:
            la, lb = lb, la
        # constant operands are folded (SAW's AIGER reader rejects AND gates with constant fan-ins):
        # AND(0,b)=0, AND(1,b)=b, XOR(0,b)=b, XOR(1,b)=~b.  Checked by validate() like everything else.
        if la <= 1:
            if op == "AND":
                r = 0 if la == 0 else lb
            else:
                r = lb ^ la
            amap[I + 1 + k] = r
            continue
        if op == "AND":
            r = mk_and(la, lb)
        elif op == "XOR":
            n1 = mk_and(la, lb)
            n2 = mk_and(la ^ 1, lb ^ 1)
            r = mk_and(n1 ^ 1, n2 ^ 1)
        else:
            raise ValueError(op)
        amap[I + 1 + k] = r
    outs = [lit(o) for o in x.outputs]
    A = len(gates)
    M = I + A
    hdr = f"aig {M} {I} 0 {len(outs)} {A}\n".encode()
    body = bytearray(hdr)
    for o in outs:
        body += f"{o}\n".encode()
    for lhs, r0, r1 in gates:
        assert lhs > r0 >= r1
        body += _enc(lhs - r0) + _enc(r0 - r1)
    return bytes(body)


def write_aig(x, path):
    with open(path, "wb") as f:
        f.write(xag_to_aig_bytes(x))


# ---------------- independent reader / simulator ----------------
def read_aig(data):
    nl = data.index(b"\n")
    hdr = data[:nl].split()
    assert hdr[0] == b"aig", hdr
    M, I, L, O, A = map(int, hdr[1:6])
    assert L == 0
    pos = nl + 1
    outs = []
    for _ in range(O):
        nl = data.index(b"\n", pos)
        outs.append(int(data[pos:nl]))
        pos = nl + 1

    def dec():
        nonlocal pos
        v, sh = 0, 0
        while True:
            c = data[pos]; pos += 1
            v |= (c & 0x7F) << sh
            sh += 7
            if not c & 0x80:
                return v
    ands = []
    for i in range(A):
        lhs = 2 * (I + L + i + 1)
        d0 = dec(); d1 = dec()
        r0 = lhs - d0; r1 = r0 - d1
        ands.append((r0, r1))
    return I, ands, outs


def simulate_aig(parsed, in_vecs, n):
    I, ands, outs = parsed
    mask = (1 << n) - 1
    val = [0] + [v & mask for v in in_vecs]
    for r0, r1 in ands:
        a = val[r0 >> 1] ^ (mask if r0 & 1 else 0)
        b = val[r1 >> 1] ^ (mask if r1 & 1 else 0)
        val.append(a & b)
    return [val[o >> 1] ^ (mask if o & 1 else 0) for o in outs]


def validate(x, aig_bytes, nvec=256, seed=1):
    """Random-vector cross-check of the AIGER encoding against verifier/sim.py (+ all-0/all-1)."""
    sys.path.insert(0, VERIFIER)
    import random
    import sim
    r = random.Random(seed)
    vecs = [r.getrandbits(nvec) for _ in range(x.n_inputs)]
    vecs = [(l & ~3) | 2 for l in vecs]          # vector 0 = all zeros, vector 1 = all ones
    a = sim.simulate(x, vecs, nvec)
    b = simulate_aig(read_aig(aig_bytes), vecs, nvec)
    return a == b
