# Interface specification and cost model (v1, 2026-11-18)

This file fixes the target functions, bit order and cost model used for every
count in the paper.

Reference implementation: `verifier/sha256_ref.py` (FIPS 180-4; constants are
derived from prime roots and checked against the published tables; the hash is
checked against hashlib and NIST vectors).

## Bit order (canonical, all interfaces)

A sequence of 32-bit words `w[0], w[1], ...` is identified with a bit vector
`v[0], v[1], ...` by

    v[i] = bit (31 - i % 32) of w[i // 32]

where bit 0 is the least significant bit. That is: words in order, each word
most-significant bit first. This equals big-endian message order: the
SHA-256 block bytes `m[0..63]` give `w[k] = m[4k]<<24 | m[4k+1]<<16 | m[4k+2]<<8 | m[4k+3]`,
so `v[i]` is bit `7 - i%8` (MSB first) of byte `m[i//8]`. The digest uses the
same rule: output bit `j` is bit `31 - j%32` of `H[j//32]`, i.e. the usual hex
digest read left to right, MSB first.

In code, `sha256_ref.words_to_bits` / `bits_to_words` implement this, and a
"packed SPEC int" (used by `sim.py`) has SPEC bit `i` at integer bit `i`.

## PRIMARY target: `primary` (512 -> 256)

    out = compress(IV, W)

* inputs `v[0..511]`: one **already-padded** 512-bit message block `W[0..15]`
  in the bit order above. The circuit does not pad; whatever block it receives
  is compressed. (For a one-block message the caller supplies the padded block;
  e.g. "abc" -> `61626380 00000000 ... 00000018`.)
* chaining value fixed to the standard IV
  `6a09e667 bb67ae85 3c6ef372 a54ff53a 510e527f 9b05688c 1f83d9ab 5be0cd19`.
* 64 standard rounds with the standard K constants, message schedule, and the
  final feed-forward `H_out[k] = H_in[k] + state[k] mod 2^32`.
* outputs `y[0..255]`: `H_out[0..7]` in the bit order above.
* For a padded single block of a message of < 56 bytes this equals the
  SHA-256 digest of that message (tested).

## SECONDARY target: `secondary` (768 -> 256)

    out = compress(CV, W)

* inputs `v[0..255]` = chaining value `CV[0..7]`; `v[256..767]` = block
  `W[0..15]`; both in the bit order above (CV first, then block).
* outputs as PRIMARY.

## Other functions recognised (not targets)

* `sha256_msg256` (256 -> 256): SHA-256 of a 32-byte message, i.e.
  `compress(IV, pad(m))` with the padding block constants baked in. Only 256
  block bits are free; NOT comparable to PRIMARY.

## Cost model (scoring)

* Basis: 2-input AND and 2-input XOR nodes; NOT is a free complemented edge;
  constants and wires free; unrestricted fanout.
* Score = AND + XOR **after** `verifier/cost.py canonicalize` (constant
  propagation, trivial identities, XOR complement normalisation, structural
  hashing, dead-node removal; no other rewriting). Always report AND and XOR
  separately, logical depth (AND and XOR each count 1), AND depth, max fanout,
  serialized size, and the canonical sha256.
* Raw counts by original gate type (incl. INV/NOT/XNOR) are reported
  separately and never mixed with canonical counts. Published counts that
  include NOT gates (e.g. NIST "Tally") are a different convention.

## Circuit interchange

* Native: `.xag` text (`verifier/xag.py write_xag/read_xag`), inputs and outputs
  in SPEC order. The canonical hash is sha256 of `write_xag(canonicalize(x))`.
* Accepted inputs: old Bristol/SCAPI, Bristol Fashion (incl. INV/EQ/EQW/MAND),
  NIST `.slp`. A non-SPEC-ordered circuit must come with an explicit interface
  map (`in_map`, `out_map`, see `sim.py`), as recorded in each
  `baselines/*/MANIFEST.json`.

## Verification status vocabulary

* `equivalent`: SAT miter UNSAT **and** the text-LRAT proof accepted by BOTH
  `verifier/lrat_check.py` (strict RUP-only checker) and the vendored
  `lrat-check`. Trust boundary: miter/Tseitin code in `equiv.py`, the XAG
  translation in `xag.py`, and `lrat_check.py`; the SAT solver is not trusted.
  The vendored `lrat-check` (drat-trim commit 2e3b2dc) is a second opinion
  only; its exit code is never used.
* `unsat_unchecked`: UNSAT without an accepted proof.
* `not_equivalent`: SAT, counterexample confirmed by simulation.
* `unverified`: timeout or error. Never a pass.
* Random-vector simulation (`sim.verify`) is a **smoke test only**.
