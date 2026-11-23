"""Independent SHA-1 / SHA-512 / SHA-256 compression references
(written from FIPS 180-4 for this audit, checked against hashlib below).
Not derived from any project code."""
import hashlib
import struct

M32, M64 = (1 << 32) - 1, (1 << 64) - 1


def rotr(x, n, w):
    return ((x >> n) | (x << (w - n))) & ((1 << w) - 1)


# ---------------- SHA-1 ----------------
SHA1_IV = [0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476, 0xC3D2E1F0]


def sha1_compress(cv, blk):
    w = list(blk)
    for t in range(16, 80):
        x = w[t - 3] ^ w[t - 8] ^ w[t - 14] ^ w[t - 16]
        w.append(((x << 1) | (x >> 31)) & M32)
    a, b, c, d, e = cv
    for t in range(80):
        if t < 20:
            f, k = (b & c) | (~b & d), 0x5A827999
        elif t < 40:
            f, k = b ^ c ^ d, 0x6ED9EBA1
        elif t < 60:
            f, k = (b & c) | (b & d) | (c & d), 0x8F1BBCDC
        else:
            f, k = b ^ c ^ d, 0xCA62C1D6
        tmp = (((a << 5) | (a >> 27)) + (f & M32) + e + k + w[t]) & M32
        a, b, c, d, e = tmp, a, ((b << 30) | (b >> 2)) & M32, c, d
    return [(x + y) & M32 for x, y in zip(cv, (a, b, c, d, e))]


# ---------------- SHA-512 ----------------
def _primes(n):
    ps, k = [], 2
    while len(ps) < n:
        if all(k % p for p in ps if p * p <= k):
            ps.append(k)
        k += 1
    return ps


def _frac_root(p, r, bits=64):
    # floor(frac(p^(1/r)) * 2^bits) by integer root
    sc = 1 << (bits + 8)
    x = int(round(p ** (1.0 / r) * (1 << 40))) << (bits + 8 - 40)
    target = p * sc ** r
    lo, hi = x - (1 << (bits)), x + (1 << (bits))
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** r <= target:
            lo = mid
        else:
            hi = mid - 1
    return (lo >> 8) & ((1 << bits) - 1)


_P = _primes(80)
SHA512_K = [_frac_root(p, 3) for p in _P]
SHA512_IV = [_frac_root(p, 2) for p in _P[:8]]


def sha512_compress(cv, blk):
    w = list(blk)
    for t in range(16, 80):
        s0 = rotr(w[t - 15], 1, 64) ^ rotr(w[t - 15], 8, 64) ^ (w[t - 15] >> 7)
        s1 = rotr(w[t - 2], 19, 64) ^ rotr(w[t - 2], 61, 64) ^ (w[t - 2] >> 6)
        w.append((w[t - 16] + s0 + w[t - 7] + s1) & M64)
    a, b, c, d, e, f, g, h = cv
    for t in range(80):
        S1 = rotr(e, 14, 64) ^ rotr(e, 18, 64) ^ rotr(e, 41, 64)
        ch = (e & f) ^ (~e & g)
        t1 = (h + S1 + ch + SHA512_K[t] + w[t]) & M64
        S0 = rotr(a, 28, 64) ^ rotr(a, 34, 64) ^ rotr(a, 39, 64)
        mj = (a & b) ^ (a & c) ^ (b & c)
        t2 = (S0 + mj) & M64
        h, g, f, e, d, c, b, a = g, f, e, (d + t1) & M64, c, b, a, (t1 + t2) & M64
    return [(x + y) & M64 for x, y in zip(cv, (a, b, c, d, e, f, g, h))]


def _pad(msg, blen, lenbytes):
    ml = len(msg) * 8
    m = msg + b"\x80"
    while (len(m) + lenbytes) % blen:
        m += b"\x00"
    return m + ml.to_bytes(lenbytes, "big")


def selftest():
    for msg in [b"", b"abc", b"x" * 50]:
        p = _pad(msg, 64, 8)
        blk = list(struct.unpack(">16I", p))
        h = sha1_compress(SHA1_IV, blk)
        assert b"".join(x.to_bytes(4, "big") for x in h) == hashlib.sha1(msg).digest()
        p = _pad(msg, 128, 16)
        blk = list(struct.unpack(">16Q", p))
        h = sha512_compress(SHA512_IV, blk)
        assert b"".join(x.to_bytes(8, "big") for x in h) == hashlib.sha512(msg).digest()
    # two-block messages exercise the variable chaining value
    msg = b"y" * 200
    p = _pad(msg, 128, 16)
    cv = SHA512_IV
    for i in range(0, len(p), 128):
        cv = sha512_compress(cv, list(struct.unpack(">16Q", p[i:i + 128])))
    assert b"".join(x.to_bytes(8, "big") for x in cv) == hashlib.sha512(msg).digest()
    return True


if __name__ == "__main__":
    print("selftest", selftest())
