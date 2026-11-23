"""port generator independent SHA-2 compression reference (FIPS 180-4), pure Python, no hashlib in the function.
Parametrised by word size: SHA-256 (32-bit, 64 rounds) and SHA-512 (64-bit, 80 rounds).
Constants are derived from prime roots with exact integer arithmetic (FIPS 180-4 4.2.2/4.2.3, 5.3.3/5.3.5);
the self-test compares with the published first/last values and with hashlib on padded messages.

SPEC bit order (same rule as verifier/SPEC.md, generalised): word sequence w[0..] -> v[i] = bit (WS-1 - i % WS)
of w[i // WS]  (words in order, each word MSB first = big-endian message order).
Inputs of the compression circuit: ALREADY-PADDED block (16 words; 512 bits for SHA-256, 1024 bits for
SHA-512), preceded, for the variable-CV targets, by 8 chaining-value words.  Output: 8 words
H_out = CV + state after all rounds (feed-forward), MSB first.  The circuit does no padding.
"""
import hashlib


def _primes(n):
    out, k = [], 2
    while len(out) < n:
        if all(k % p for p in out if p * p <= k):
            out.append(k)
        k += 1
    return out


def _iroot(x, r):
    lo, hi = 0, 1
    while hi ** r <= x:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** r <= x:
            lo = mid
        else:
            hi = mid
    return lo


def _frac(p, r, ws):
    return _iroot(p << (ws * r), r) & ((1 << ws) - 1)


class Sha2:
    def __init__(s, ws):
        assert ws in (32, 64)
        s.ws = ws; s.M = (1 << ws) - 1
        s.R = 64 if ws == 32 else 80
        s.IV = tuple(_frac(p, 2, ws) for p in _primes(8))
        s.K = tuple(_frac(p, 3, ws) for p in _primes(s.R))
        if ws == 32:
            s.bs0, s.bs1, s.ss0, s.ss1 = (2, 13, 22), (6, 11, 25), (7, 18, 3), (17, 19, 10)
        else:
            s.bs0, s.bs1, s.ss0, s.ss1 = (28, 34, 39), (14, 18, 41), (1, 8, 7), (19, 61, 6)

    def rotr(s, x, n): return ((x >> n) | (x << (s.ws - n))) & s.M
    def Bsig0(s, x): a, b, c = s.bs0; return s.rotr(x, a) ^ s.rotr(x, b) ^ s.rotr(x, c)
    def Bsig1(s, x): a, b, c = s.bs1; return s.rotr(x, a) ^ s.rotr(x, b) ^ s.rotr(x, c)
    def ssig0(s, x): a, b, c = s.ss0; return s.rotr(x, a) ^ s.rotr(x, b) ^ (x >> c)
    def ssig1(s, x): a, b, c = s.ss1; return s.rotr(x, a) ^ s.rotr(x, b) ^ (x >> c)
    def Ch(s, e, f, g): return (e & f) ^ (~e & s.M & g)
    def Maj(s, a, b, c): return (a & b) ^ (a & c) ^ (b & c)

    def schedule(s, W16):
        W = list(W16)
        for t in range(16, s.R):
            W.append((s.ssig1(W[t-2]) + W[t-7] + s.ssig0(W[t-15]) + W[t-16]) & s.M)
        return W

    def compress(s, H, W16):
        W = s.schedule(W16)
        a, b, c, d, e, f, g, h = H
        for t in range(s.R):
            t1 = (h + s.Bsig1(e) + s.Ch(e, f, g) + s.K[t] + W[t]) & s.M
            t2 = (s.Bsig0(a) + s.Maj(a, b, c)) & s.M
            h, g, f, e, d, c, b, a = g, f, e, (d + t1) & s.M, c, b, a, (t1 + t2) & s.M
        return [(x + y) & s.M for x, y in zip(H, (a, b, c, d, e, f, g, h))]

    def pad(s, msg):
        bs = 16 * s.ws // 8
        lb = 8 if s.ws == 32 else 16
        m = msg + b"\x80"
        m += b"\x00" * ((-len(m) - lb) % bs)
        m += (8 * len(msg)).to_bytes(lb, "big")
        wb = s.ws // 8
        return [[int.from_bytes(m[i + wb * k:i + wb * (k + 1)], "big") for k in range(16)] for i in range(0, len(m), bs)]

    def hash(s, msg):
        H = list(s.IV)
        for blk in s.pad(msg):
            H = s.compress(H, blk)
        return b"".join(w.to_bytes(s.ws // 8, "big") for w in H)

    def words_to_bits(s, words):
        return [(w >> (s.ws - 1 - j)) & 1 for w in words for j in range(s.ws)]

    def bits_to_words(s, bits):
        return [sum(bits[s.ws * k + j] << (s.ws - 1 - j) for j in range(s.ws)) for k in range(len(bits) // s.ws)]

    def words_to_int(s, words):
        """packed SPEC int: SPEC bit i at integer bit i"""
        v = 0
        for i, b in enumerate(s.words_to_bits(words)):
            v |= b << i
        return v


SHA256, SHA512 = Sha2(32), Sha2(64)

if __name__ == "__main__":
    import os
    assert SHA512.IV[0] == 0x6a09e667f3bcc908 and SHA512.IV[7] == 0x5be0cd19137e2179
    assert SHA512.K[0] == 0x428a2f98d728ae22 and SHA512.K[79] == 0x6c44198c4a475817
    assert SHA256.IV[0] == 0x6a09e667 and SHA256.K[63] == 0xc67178f2
    for S, hf in ((SHA256, hashlib.sha256), (SHA512, hashlib.sha512)):
        for n in list(range(0, 300)) + [1000, 4096]:
            m = os.urandom(n)
            assert S.hash(m) == hf(m).digest(), (S.ws, n)
    print("sha2ref: constants + hash vs hashlib OK")
