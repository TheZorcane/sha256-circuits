"""Independent SHA-256 reference (FIPS 180-4), pure Python, no hashlib.

Constants are *derived* (fractional parts of square / cube roots of the first
primes, FIPS 180-4 sections 4.2.2 and 5.3.3) using exact integer arithmetic, so
no table can be mistyped; tests compare them with the published values and the
whole function with hashlib and NIST test vectors.

Public API
  compress(H_in, W_block)   -> 8 x u32   one compression incl. feed-forward
  compress_fixed_iv(W_block) -> 8 x u32  compress(IV, W_block)   (SPEC PRIMARY)
  sha256(msg: bytes) -> bytes            full hash with padding (for tests)
  pad(msg: bytes) -> list[list[u32]]     padded message as 16-word blocks
"""

MASK = 0xFFFFFFFF


def _primes(n):
    out, k = [], 2
    while len(out) < n:
        if all(k % p for p in out if p * p <= k):
            out.append(k)
        k += 1
    return out


def _iroot(x, r):
    """floor(x ** (1/r)) for non-negative integer x, exact."""
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


def _frac32(p, r):
    """First 32 bits of the fractional part of p ** (1/r)."""
    return _iroot(p << (32 * r), r) & MASK


IV = tuple(_frac32(p, 2) for p in _primes(8))
K = tuple(_frac32(p, 3) for p in _primes(64))


def _rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK


def _bsig0(x):
    return _rotr(x, 2) ^ _rotr(x, 13) ^ _rotr(x, 22)


def _bsig1(x):
    return _rotr(x, 6) ^ _rotr(x, 11) ^ _rotr(x, 25)


def _ssig0(x):
    return _rotr(x, 7) ^ _rotr(x, 18) ^ (x >> 3)


def _ssig1(x):
    return _rotr(x, 17) ^ _rotr(x, 19) ^ (x >> 10)


def _ch(x, y, z):
    return (x & y) ^ (~x & z & MASK)


def _maj(x, y, z):
    return (x & y) ^ (x & z) ^ (y & z)


def compress(H_in, W_block):
    """SHA-256 compression function: 64 rounds + feed-forward (FIPS 180-4 6.2.2)."""
    H_in = [int(h) for h in H_in]
    W_block = [int(w) for w in W_block]
    if len(H_in) != 8 or len(W_block) != 16:
        raise ValueError("need 8 chaining words and 16 block words")
    if any(not 0 <= v <= MASK for v in H_in + W_block):
        raise ValueError("words must be u32")
    W = list(W_block)
    for t in range(16, 64):
        W.append((_ssig1(W[t - 2]) + W[t - 7] + _ssig0(W[t - 15]) + W[t - 16]) & MASK)
    a, b, c, d, e, f, g, h = H_in
    for t in range(64):
        t1 = (h + _bsig1(e) + _ch(e, f, g) + K[t] + W[t]) & MASK
        t2 = (_bsig0(a) + _maj(a, b, c)) & MASK
        h, g, f, e, d, c, b, a = g, f, e, (d + t1) & MASK, c, b, a, (t1 + t2) & MASK
    return [(x + y) & MASK for x, y in zip(H_in, (a, b, c, d, e, f, g, h))]


def compress_fixed_iv(W_block):
    return compress(IV, W_block)


def pad(msg: bytes):
    """FIPS 180-4 5.1.1 padding; returns list of 16-word (big-endian) blocks."""
    ml = len(msg) * 8
    m = bytes(msg) + b"\x80"
    m += b"\x00" * ((56 - len(m)) % 64)
    m += ml.to_bytes(8, "big")
    words = [int.from_bytes(m[i:i + 4], "big") for i in range(0, len(m), 4)]
    return [words[i:i + 16] for i in range(0, len(words), 16)]


def sha256(msg: bytes) -> bytes:
    H = list(IV)
    for blk in pad(msg):
        H = compress(H, blk)
    return b"".join(h.to_bytes(4, "big") for h in H)


# ---- SPEC bit order helpers (see verifier/SPEC.md) ----------------------------
# A sequence of u32 words <-> a bit vector: bit i = bit (31 - i % 32) of word i // 32.

def words_to_bits(words):
    return [(words[i // 32] >> (31 - i % 32)) & 1 for i in range(32 * len(words))]


def bits_to_words(bits):
    if len(bits) % 32:
        raise ValueError("bit length must be a multiple of 32")
    out = []
    for w in range(len(bits) // 32):
        v = 0
        for j in range(32):
            v = (v << 1) | (bits[32 * w + j] & 1)
        out.append(v)
    return out


def words_to_int(words):
    """SPEC bit vector packed as an integer whose bit i is SPEC bit i."""
    v = 0
    for i, b in enumerate(words_to_bits(words)):
        v |= b << i
    return v


def primary(bits512):
    """SPEC PRIMARY on a 512-entry bit list -> 256-entry bit list."""
    return words_to_bits(compress_fixed_iv(bits_to_words(bits512)))


def secondary(bits768):
    """SPEC SECONDARY: bits 0..255 = chaining value, 256..767 = block."""
    w = bits_to_words(bits768)
    return words_to_bits(compress(w[:8], w[8:]))
