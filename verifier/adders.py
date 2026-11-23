"""Small reference circuits for harness demos/tests: n-bit adders mod 2^n.

Inputs: a bits 0..n-1 (LSB first), b bits 0..n-1; outputs: sum bits LSB first.
"""

from xag import XAG


def ripple_maj1(n):
    """Ripple carry; carry = ((a^c)&(b^c))^c  (1 AND per bit)."""
    x = XAG(2 * n)
    c, out = 0, []
    for i in range(n):
        a, b = x.input(i), x.input(n + i)
        axb = x.XOR(a, b)
        out.append(x.XOR(axb, c))
        if i < n - 1:
            c = x.XOR(x.AND(x.XOR(a, c), x.XOR(b, c)), c)
    x.outputs = out
    return x


def ripple_maj2(n):
    """Ripple carry; carry = (a&b) ^ (c&(a^b))  (2 ANDs per bit)."""
    x = XAG(2 * n)
    c, out = 0, []
    for i in range(n):
        a, b = x.input(i), x.input(n + i)
        p = x.XOR(a, b)
        out.append(x.XOR(p, c))
        if i < n - 1:
            c = x.XOR(x.AND(a, b), x.AND(c, p))
    x.outputs = out
    return x


def kogge_stone(n):
    """Parallel-prefix (Kogge-Stone) adder."""
    x = XAG(2 * n)
    a = [x.input(i) for i in range(n)]
    b = [x.input(n + i) for i in range(n)]
    p = [x.XOR(a[i], b[i]) for i in range(n)]
    g = [x.AND(a[i], b[i]) for i in range(n)]
    G, P = g[:], p[:]
    d = 1
    while d < n:
        G2, P2 = G[:], P[:]
        for i in range(d, n):
            # (G,P)_i o (G,P)_{i-d}: G = G_i | (P_i & G_{i-d}); G_i, P_i&G_{i-d}
            # are disjoint, so OR == XOR
            G2[i] = x.XOR(G[i], x.AND(P[i], G[i - d]))
            P2[i] = x.AND(P[i], P[i - d])
        G, P = G2, P2
        d *= 2
    x.outputs = [p[0]] + [x.XOR(p[i], G[i - 1]) for i in range(1, n)]
    return x


def broken(n, bit=None):
    """ripple_maj2 with the AND in one carry replaced by XOR (a real bug)."""
    bit = n // 2 if bit is None else bit
    x = XAG(2 * n)
    c, out = 0, []
    for i in range(n):
        a, b = x.input(i), x.input(n + i)
        p = x.XOR(a, b)
        out.append(x.XOR(p, c))
        if i < n - 1:
            ab = x.XOR(a, b) if i == bit else x.AND(a, b)
            c = x.XOR(ab, x.AND(c, p))
    x.outputs = out
    return x
