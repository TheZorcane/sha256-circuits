"""port generator independent SHA-1 compression reference (FIPS 180-4 section 6.1), pure Python.
SPEC (frozen here): input = ONE already-padded 512-bit block, 16 words big-endian, word-major MSB-first bit order
(v[i] = bit (31 - i%32) of W[i//32]); the chaining value is the fixed standard IV
(67452301 efcdab89 98badcfe 10325476 c3d2e1f0); output = H_out[0..4] = IV + state after 80 rounds, same bit order
(160 bits).  The circuit does not pad.  Variable-CV variant (not scored): 160 CV bits first, then the block."""
import hashlib
M = 0xFFFFFFFF
IV = (0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476, 0xC3D2E1F0)
K = (0x5A827999, 0x6ED9EBA1, 0x8F1BBCDC, 0xCA62C1D6)


def rotl(x, n): return ((x << n) | (x >> (32 - n))) & M


def f(t, b, c, d):
    if t < 20: return (b & c) | (~b & M & d)
    if t < 40: return b ^ c ^ d
    if t < 60: return (b & c) | (b & d) | (c & d)
    return b ^ c ^ d


def schedule(W16):
    W = list(W16)
    for t in range(16, 80):
        W.append(rotl(W[t-3] ^ W[t-8] ^ W[t-14] ^ W[t-16], 1))
    return W


def compress(H, W16):
    W = schedule(W16)
    a, b, c, d, e = H
    for t in range(80):
        T = (rotl(a, 5) + f(t, b, c, d) + e + K[t // 20] + W[t]) & M
        e, d, c, b, a = d, c, rotl(b, 30), a, T
    return [(x + y) & M for x, y in zip(H, (a, b, c, d, e))]


def pad(msg):
    m = msg + b"\x80"
    m += b"\x00" * ((-len(m) - 8) % 64)
    m += (8 * len(msg)).to_bytes(8, "big")
    return [[int.from_bytes(m[i + 4 * k:i + 4 * k + 4], "big") for k in range(16)] for i in range(0, len(m), 64)]


def sha1(msg):
    H = list(IV)
    for b in pad(msg): H = compress(H, b)
    return b"".join(w.to_bytes(4, "big") for w in H)


def words_to_bits(words): return [(w >> (31 - j)) & 1 for w in words for j in range(32)]
def bits_to_words(bits): return [sum(bits[32 * k + j] << (31 - j) for j in range(32)) for k in range(len(bits) // 32)]
def words_to_int(words):
    v = 0
    for i, b in enumerate(words_to_bits(words)): v |= b << i
    return v


if __name__ == "__main__":
    import os
    for n in list(range(0, 200)) + [1000]:
        m = os.urandom(n)
        assert sha1(m) == hashlib.sha1(m).digest(), n
    print("sha1ref OK vs hashlib")
