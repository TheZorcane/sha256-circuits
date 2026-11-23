# Artifact: Smaller verified Boolean circuits for SHA-256

Ido Rosenbaum, 2026. License: see `LICENSE` (CC BY 4.0 for documentation, circuits and data; MIT for the code;
the vendored drat-trim sources in `verifier/vendor/` keep their own license, `verifier/vendor/LICENSE`).

- Paper: Ido Rosenbaum, *Smaller verified Boolean circuits for SHA-256*, IACR Cryptology ePrint Archive
  2026/EPRINT_ID.
- Archived version of this artifact: https://doi.org/ZENODO_DOI
- Cite: see `CITATION.cff`.

Quick start: `sh reproduce.sh` (about 7 minutes on a laptop; needs Python 3, a C compiler and CaDiCaL 3.x;
compare `out/summary.txt` with `expected/summary.txt`).

This artifact contains the circuits of the paper, the verifier, the generators, the proof manifests, the SAW
and Cryptol scripts that link the verifier's reference to Galois's Cryptol FIPS 180-4 specification, and
scripts that reproduce the results.

## Circuits (fixed-IV SHA-256 compression; canonical AND + XOR, NOT free)

| file | AND | XOR | total |
|---|---|---|---|
| `circuits/sha256_compress_96273.xag` (best total) | 21,783 | 74,490 | 96,273 |
| `circuits/sha256_compress_96324.xag` (trade-off) | 21,746 | 74,578 | 96,324 |
| `circuits/sha256_compress_100658.xag` (fewest AND of ours) | 21,745 | 78,913 | 100,658 |

`circuits/bristol/` has the same three circuits in Bristol Fashion format; `circuits/proofbase/` has the two
proof bases (99,630 and 99,669 gates) proved directly against the reference, and the first linear-layer pass
(96,274); `circuits/ports/` has the SHA-1 circuit (49,851 gates, 10,552 AND) and the variable-CV SHA-256 circuit
(103,669 gates, 22,293 AND; 768 inputs: the 8 chaining-value words, then the 16 block words).

Interface: inputs in the order of `verifier/SPEC.md` (512 block bits, words W_0..W_15, each MSB first), outputs
the 256 bits of H_0..H_7, MSB first. The Bristol Fashion files use one input group of 512 wires and one output
group of 256 wires in the same order (the Bristol Fashion collection's own SHA-256 uses a different word and bit
order).

## What is where

| path | content |
|---|---|
| `reproduce.sh` | quick mode (default) and `full` mode; PASS/FAIL per step in `out/summary.txt` |
| `expected/` | expected `summary.txt` and logs of the quick run |
| `SHA256SUMS` | hashes of every circuit file and generator data file |
| `verifier/` | XAG library and parsers (`xag.py`), canonical cost model (`cost.py`, `SPEC.md`), simulator (`sim.py`), FIPS 180-4 reference (`sha256_ref.py`), SAT-miter harness (`equiv.py`), strict RUP-only LRAT checker (`lrat_check.py`), vendored drat-trim `lrat-check` source (secondary check; no claim depends on it) |
| `compositional_proof/` | lemma prover (`lemmas.py`), cut finder (`cuts.py`, `cutmap.py`), word-level reference (`sharef.py`, with the single-tail and split-tail variants), composition checker (`check_compose.py`), run merger (`merge_runs.py`), simulation and mutation harness (`attack.py`; not used by `reproduce.sh`) |
| `twin_cut_certificate/` | LRAT certificate for linear-layer rewrites (`run_certificate.sh` wraps `prep.py`, `run.py`, `compose.py`; `negctl.py` negative controls; `check3.py`, a third plane-matching checker) |
| `rewrite_checks/` | GF(2) checker `lincheck.py`, independent re-implementation `independent/check2.py`, structural `rewrite_check.py`, negative controls `mutants.py` |
| `fips_spec_proof/` | SAW 1.6 / Cryptol proofs that the verifier's reference equals Galois's Cryptol SHA-256 specification, and `check_ref_link.py` (see below) |
| `generator/` | circuit generators: `generate_sha256.py` (96,273 / 96,324 lines), `generate_sha256_joint.py` (100,658), shared modules, `peep.py`, the joint sigma program `joint_sigma_slp_105.txt`; SHA-1/SHA-2 ports `gen_sha1.py`, `gen_sha2.py` with references and simulators |
| `residual_walk/` | deterministic linear-layer search (`walk.sh`, C sources) |
| `comparators/` | textbook generator (`textbook.py`, with helper module `rescore.py`), emp-toolkit transliteration (`emp_sha256.py`), independent SHA-1/SHA-512 references (`refs.py`), public-file scorer and download script (`score_public.py`, `fetch_public_baselines.sh`), `RAW_SHA256SUMS` of the scored public files, ABC control recipe (`abc_control/`: `abc_control.py`, `xag.genlib`, BLIF converters `xagio.py`) |
| `manifests/` | proof manifests: per-lemma JSON (CNF hash, proof hash, solver, checker verdicts, timings) and composition reports for both proof bases and for the 100,658 circuit; twin-cut certificate logs and per-batch manifests; the two tail-lemma records used by quick mode |
| `tools/` | `check_circuits.py` (recount incl. the NOT-counting tally, simulate, export and re-read Bristol Fashion), `negative_controls.sh`, `mutate_and.py` (negative control for the compositional prover) |

## Dependencies
- Python 3 (standard library only; tested with 3.12 and 3.15).
- A C compiler (`cc`).
- CaDiCaL 3.x on `PATH` as `cadical` (tested with 3.1.0). CaDiCaL 2.2.1 occasionally emits LRAT steps with a
  redundant hint, which the strict checker rejects (fail-safe); use 3.x.
- For `fips_spec_proof/` only: SAW 1.6 (with Cryptol 3.6), yices 2.8 on `PATH`, and a checkout of
  https://github.com/GaloisInc/cryptol-specs at commit 8b74f5fc2b77fc7bbb587ba2217990b1963e9d6c.

## Reproduction

```sh
sh reproduce.sh                 # quick mode: about 7 minutes, up to 4 threads, < 2 GB disk
JOBS=2 sh reproduce.sh          # fewer SAT jobs
sh reproduce.sh full            # quick mode + rebuilds (about 30 min) + all compositional proofs
                                # (about 1 hour and 20 GB of free disk per proof base)
# link from the reference to the Cryptol FIPS 180-4 specification
# (about 10 minutes with 3 SAW processes, the default; about 25 minutes with SAW_JOBS=1):
CRYPTOL_SPECS=/path/to/cryptol-specs SAW=/path/to/saw bash fips_spec_proof/run.sh
```
`reproduce.sh` and `fips_spec_proof/run.sh` exit with status 0 only if every required check passed; compare
`out/summary.txt` with `expected/summary.txt`, and `fips_spec_proof/out/summary.txt` with
`fips_spec_proof/expected/summary.txt` (all ten SAW logs, `gen.log`, `env.txt` and `ref_link.log` of our run are
in `fips_spec_proof/expected/`; paths in them are shortened to `<artifact>`). `python3
fips_spec_proof/summarize.py fips_spec_proof/expected` re-checks the shipped logs without SAW.

### What quick mode does
1. `sha256sums`, `counts` (canonical counts, depth and the NOT-counting tally 99,035 / 99,054),
   `simulation` (100,000 random blocks per circuit plus 4 fixed inputs: for the fixed-IV circuits the padded
   empty and "abc" messages and the all-zero and all-one blocks; the variable-CV circuit is included),
   `bristol_reread`.
2. `gf2_checks`: the two independent GF(2) checkers on (proof base, final circuit) for both lines; all four
   checks must print PASS. `lincheck.py` also counts the matched AND gates whose affine correction is zero, i.e.
   that compute exactly the function of their proof-base AND (21,783 of 21,783 and 21,746 of 21,746); output
   equality does not depend on this count.
3. `twin_cut_96273`, `twin_cut_96324`: LRAT certificates final circuit = proof base (about 40 s each).
4. `compositional_96273_108of110`: re-proves all 46 schedule and 62 round lemmas of the 96,273 proof base
   against the reference (about 2 minutes with 4 jobs), merges them with the recorded manifests of the two
   fused-round tail lemmas, and runs the composition checker, which rebuilds those two cones and checks their
   hashes against the manifests. **Only this proof base is re-proved in quick mode**; for the 96,324 proof base
   and the 100,658 circuit quick mode relies on the recorded manifests in `manifests/` (full mode re-proves
   them).
   `compositional_negative_control`: `tools/mutate_and.py` complements one input of one AND gate in the
   round-10 cone of the 96,273 proof base; `lemmas.py` must report the round-10 lemma `not_equivalent`
   (a word bit of the reference is then computed by no node, and a whole-circuit simulation against the
   reference confirms a counterexample).
5. `ref_link_329`: `fips_spec_proof/check_ref_link.py all` (see below), without SAW.
6. `negative_controls`: `tools/negative_controls.sh` (40 mutated twin-cut lemmas must be refuted; 21 mutants
   must be rejected by both GF(2) checkers).
7. `variable_cv`: regenerates the variable-CV generator output (22,293 AND / 104,023 gates), simulates it, and
   checks the shipped 103,669-gate circuit against it with both GF(2) checkers (it differs only in the linear
   layer: the joint sigma pass). The step fails if the simulator (`generator/simsha2.py`, which exits nonzero on
   any mismatch) or either GF(2) checker fails; each exit status is checked separately.
8. `textbook_baseline` (fails if any of its simulations fails), `emp_transliteration`.

The simulators `generator/simsha1.py`, `generator/simsha2.py` and `generator/simcheck.py` exit nonzero on any
mismatch, so the full-mode `ports_sha1_sha512` step also fails on a simulation mismatch.

### Public comparator files (manual fetch, needs network)
The public circuits are not shipped. To score them:
```sh
sh comparators/fetch_public_baselines.sh public_raw    # curl from usnistgov/Circuits (NIST, both files),
                                                       # matteocam/sha256-circuit (Campanelli et al.) and
                                                       # nigelsmart.github.io/MPC-Circuits (Bristol Fashion)
python3 comparators/score_public.py public_raw         # or: PUBLIC_RAW=public_raw sh reproduce.sh
python3 comparators/abc_control/abc_control.py public_raw --abc /path/to/abc [--deepsyn]   # ABC controls
```
`score_public.py` checks each file against `comparators/RAW_SHA256SUMS`, prints its own tally, the direct and the
canonical count (NIST variable-CV: file tally 118,287, canonical 116,239 = 22,573 AND + 93,666 XOR), and for the
two variable-CV files also the IV-folded count (NIST 113,386; Bristol Fashion 129,422). If a file has changed
upstream its hash no longer matches and it is reported, not scored.

ABC control comparators (paper, Table 2, rows labelled "ABC control"). `abc_control.py` puts the NIST fixed-IV
and Campanelli et al. files in SPEC order, canonicalises them, writes BLIF, runs ABC (`strash`, then nothing,
the `compress2rs` script once, or `&deepsyn -T 600 -S 1`, then `read_genlib xag.genlib; map -a`), reads the
mapped BLIF back, canonicalises, scores and simulates it, and compares with the counts and canonical SHA-256
of our historical run (2026-11-18; ABC commit c2c7046d built with the default Makefile and
`ABC_USE_NO_READLINE=1` in an Ubuntu 24.04 container, gcc, on x86-64 (AMD EPYC 9B14, Google Cloud c3d); one run per
pipeline, single-threaded). Both raw files must be present; a missing or altered file, an ABC failure, a failed
simulation or an unexpected result gives exit status 1.

| input | pipeline | AND | XOR | total | canonical SHA-256 |
|---|---|---|---|---|---|
| NIST | map -a | 22,296 | 89,152 | 111,448 | `17902aee486deaf82f3827ce37e7e4aa45ed19a6c5d830f7f91095d612c2d012` |
| NIST | compress2rs; map -a | 62,580 | 48,348 | 110,928 | `380d52dba08153a00e8da7794360a9df838211061c6e461edce3ecc1e67a8665` |
| NIST | deepsyn; map -a | 61,978 | 48,336 | 110,314 | `bdba2a97d5547aa1d7324c22d1a717c2a2bb26b47a7affc850401943eaff4f8a` |
| Campanelli et al. | map -a (historical, Linux build) | 22,272 | 89,261 | 111,533 | `1fa2ec35964b19effc1802cfe5638390701ec881ff81eb2e8482e9477e840af7` |
| Campanelli et al. | map -a (macOS arm64 build, Apple clang 21) | 22,272 | 89,224 | 111,496 | `8fd4a25dac7d37858ed50b77d1d8c99c502d8e3656ba998f273f783bd74869a5` |
| Campanelli et al. | compress2rs; map -a | 61,356 | 48,424 | **109,780** | `b68730ca7fe99898c1a2a05b2bba5628a2324e9809d0dd9764a0a3287f2d9c72` |
| Campanelli et al. | deepsyn; map -a | 62,086 | 48,325 | 110,411 | `f957fb56a5334193e958f02a238d20bacf5b507602b0f5f8ff4e1854317e3041` |

Each matched the reference on 100,004 simulated inputs in our run; they are smoke-tested, not proved. ABC is not
shipped, and `reproduce.sh` does not run this step.

The Campanelli et al. map-only row depends on the ABC build. On the byte-identical input BLIF, a Linux arm64
build of c2c7046d (Ubuntu 24.04 container, gcc 13.3.0) reproduces all four map-only and compress2rs rows with
the historical canonical hashes. A macOS arm64 build of the same commit (Apple clang 21.0.0, macOS 26.7) also
reproduces the NIST rows and the headline 109,780 row, but maps Campanelli et al. with `map -a` alone to 111,496
gates (second row above); both mappings pass simulation. The cause is not established. `abc_control.py` accepts
either recorded result for that row and says which one it matched; every other row must match the historical
hash. We could not build ABC for x86-64 locally (the emulated compiler crashed), so the x86-64 result itself
was not re-run. The ABC flows trade XOR for AND (about
62,000 AND instead of about 22,300), so the paper reports AND and total for every row.

Not reproduced by the artifact: the zk.golf measurements (their sources carry no license), Testa et al.
(figures from the paper), and the intermediate circuits of the paper's Table 3.

### Full mode (additional steps)
| step | what | time / resources |
|---|---|---|
| `rebuild_generators` | regenerate both proof bases and the 100,658 circuit from scratch; must be byte-identical | 12 min wall (3 processes) |
| `rebuild_walks` | rewrite + linear-layer walks; final circuits must be byte-identical | about 18 min, 1 thread |
| `compositional_96273_full`, `compositional_96324_full` | all 110 lemmas including the two tail lemmas | about 1 hour and about 20 GB of free disk each (tail LRAT 1.9-8.6 GB, solved in 4-25 min and checked by the strict checker in 2.5-12 min) |
| `compositional_100658_full` | 109 lemmas, single-tail mode | about 15 min, about 6 GB |
| `ports_sha1_sha512` | regenerate SHA-1 (byte-identical) and SHA-512 (counts), simulate against independent references | about 4 min |

### Link to the Cryptol FIPS 180-4 specification (`fips_spec_proof/`)
`run.sh` imports the verifier's own reference generator `compositional_proof/sharef.py` for each of the three
verifier variants (`plain`: 120 lemmas; `tail`, flag `VERIFY_TAIL=1`: 109 lemmas, used for the 100,658 circuit;
`split`, flags `VERIFY_TAIL=1 VERIFY_SPLIT=1`: 110 lemmas, used for the two proof bases), writes every reference
fragment as an AIGER file (validated by simulation against `verifier/sim.py`), generates Cryptol modules for the
fragments and their composition from the verifier's lemma and wiring tables, and runs SAW (`SAW_JOBS`
processes at once, default 3):

- 678 fragment theorems (339 per style): each reference fragment equals a Cryptol function built from the
  specification's word operations (style S: operations written as in `SHA2/Specification.cry`; style I:
  `SHA2Internal`'s exported operations);
- 3 composition theorems: for every 512-bit block b, the composition equals
  `SHA2Internal::processBlock_Common` with the standard IV applied to b;
- 3 x 449 bridge theorems: for every message length L = 0..447 (447 twice), the composition applied to the
  padding of an L-bit message M as specified in FIPS 180-4 (restated in `saw/bridge_lengths.saw`, because the
  specification's `pad` is private) equals `SHA2::Specification::hash M`;
- the FIPS "abc" vector through every function, two concrete vectors per fragment, and a negative control (the
  round-10 fragment with one AND changed to XOR), which must be rejected with a genuine counterexample;
- finally `check_ref_link.py all build`: for all 329 lemmas of the three circuits (110 + 110 + 109), the reference
  side of the SAT miter built by `lemmas.py` ends with exactly the gates of the SAW fragment, node for node, under
  a substitution that maps each fragment input bit to one literal of the miter (IV-valued words are the same
  constants in both); the miter's outputs are the images of the fragment's outputs; and the `.aig` file that
  SAW read is byte-identical to the converted fragment. With `REF_LINK_SELFTEST=1` one gate of every fragment is
  flipped and every lemma must fail.

Success criterion: the exit status of `run.sh`. It is 0 only if the `cryptol-specs` checkout is at commit
8b74f5fc2b77 (checked before any proof runs), every SAW run exits 0, the negative control rejects the mutant,
`check_ref_link.py` succeeds, and `summarize.py out` succeeds. `summarize.py` requires the exact six fragment
logs with 120 / 109 / 110 distinct theorems each (the same fragment names in styles S and I), the three
composition logs each with `W_I`, `B_S`, `B_S_0` .. `B_S_447` and `B_S_ALL_LENGTHS` (451 distinct `PROVED`
lines) and `sanity OK`, the negative-control log, the 329 per-lemma reference links with the three per-circuit
lines and the final `ALL REF LINKS OK: 329 lemmas`, `gen.log` and the spec commit in `env.txt`. A missing log, a
missing or repeated theorem, a nonzero SAW exit, a timeout or a failure line gives exit status 1. We tested this
logic on the shipped logs (exit 0) and on broken copies of them (a missing fragment, composition or link log, a
failed theorem, a wrong count, a missing bridge length, a nonzero SAW exit, a timeout, a failed negative control,
a failed link line, a wrong or unknown commit, a missing `env.txt`: exit 1 in every case), and ran `run.sh` with
stand-ins for SAW (a failing one: exit 1; one that replays the shipped logs: exit 0) and with a checkout that is
not at the pinned commit (exit 1). SAW itself was not re-run for this revision of the scripts.

Within `cryptol-specs`, `SHA2Internal` is related to `SHA2::Specification` only by the tests in
`Primitive/Keyless/Hash/SHA2Internal/Equivalence.cry` (run by the repository's documentation tests): exhaustive
on all 10-bit messages and all 1-byte messages, random on 256-bit and 256-byte messages.

Files: `gen.py` (generator), `aiger.py` (XAG to AIGER converter and an independent AIGER reader used to validate
it), `cry/Frag.cry` (fragment functions; untrusted intermediates, validated by the composition and bridge
theorems), `saw/compose.saw`, `saw/bridge_lengths.saw`, `saw/negative.saw`, `run_saw.py` (runs SAW with a time
and memory budget), `summarize.py`, `check_ref_link.py`. A direct proof that the two compositions (styles S and
I) are equal on all blocks was not attempted to completion; neither chain needs it. A monolithic SAW check of the
96,273 circuit against the specification, without the decomposition, did not finish within 2 hours and is not
included.

Trust base of this step: SAW, Cryptol, What4 and yices (no solver certificates), plus `gen.py`, `aiger.py` and
`check_ref_link.py` (tested, not proved).

## Exact commands
Run from this directory with `export PYTHONPATH=$PWD/verifier`.

Generators (from `generator/`):
```sh
J=$(python3 -c "print(','.join(map(str,range(1,62))))")
# 96,273 line (proof base 99,630)
GEN_FUSEF=1 GEN_JS=joint_sigma_slp_105.txt MADD2_REG=0 GEN_LAMBDA=8 python3 generate_sha256.py a_B.xag $J
python3 peep.py a_B.xag a_Bp.xag                       # = circuits/proofbase/sha256_compress_96273_proofbase.xag
# 96,324 line (proof base 99,669)
GEN_FUSEF=1 GEN_JS=joint_sigma_slp_105.txt MADD2_REG=0 GEN_LAMBDA=8 GEN_MOVE0=1 GEN_J12=1,2 python3 generate_sha256.py b_B.xag $J
python3 peep.py b_B.xag b_Bp.xag                       # = circuits/proofbase/sha256_compress_96324_proofbase.xag
# 100,658 (fewest AND)
JOINT_LAMBDA=8 python3 generate_sha256_joint.py c_raw.xag final+wfuse+sig105
python3 peep.py c_raw.xag c.xag                        # = circuits/sha256_compress_100658.xag
# variable-CV SHA-256 generator output (22,293 AND / 104,023 gates)
python3 gen_sha2.py v.xag 32 0
```
Switches: `GEN_LAMBDA` (weight of ANDs in the planner objective), `GEN_FUSEF` (difference-word chains in the
fused rounds 62/63), `GEN_DELTA` (self-delta by-product, default on), `GEN_T1F` (chains in T_1, default on),
`GEN_MOVE0`, `GEN_J12` (joint column 0 / joint rounds), `MADD2_REG=0` (proof variant: consumers rebuild
difference words with their own XORs, so that word cuts exist for the compositional proof). Planner caches go
to `GEN_CACHE` (default `generator/cache`; created if missing).

Linear-layer passes:
```sh
python3 rewrite_checks/rewrite_check.py a_Bp.xag - a_Bpr.xag
sh residual_walk/walk.sh a_Bpr.xag a_Bp.xag a_pass1.xag                       # 96,274
sh residual_walk/walk.sh a_pass1.xag a_Bp.xag a_pass2.xag 100000000 8          # 96,273
python3 rewrite_checks/rewrite_check.py b_Bp.xag - b_Bpr.xag
sh residual_walk/walk.sh b_Bpr.xag b_Bp.xag b_pass1.xag                       # 96,324
```
(`rewrite_check.py` with `-` as candidate prints `DIFFERENT (no claim)` and writes the rewritten circuit.)

Proofs:
```sh
# compositional proof of a proof base (split-tail mode; use VERIFY_SPLIT=0 for the 100,658 circuit)
cd compositional_proof
VERIFY_TAIL=1 VERIFY_SPLIT=1 python3 lemmas.py ../work/pb.xag ../work/L --jobs 4 --timeout 10800
VERIFY_TAIL=1 VERIFY_SPLIT=1 python3 check_compose.py ../work/pb.xag ../work/L/lemmas.json ../work/L
cd ..
# twin-cut certificate final == proof base
sh twin_cut_certificate/run_certificate.sh circuits/proofbase/sha256_compress_96273_proofbase.xag circuits/sha256_compress_96273.xag out/tc 4
# algebraic checks
python3 rewrite_checks/lincheck.py circuits/proofbase/sha256_compress_96273_proofbase.xag circuits/sha256_compress_96273.xag
(cd rewrite_checks/independent && python3 check2.py ../../circuits/proofbase/sha256_compress_96273_proofbase.xag ../../circuits/sha256_compress_96273.xag)
# negative controls
sh tools/negative_controls.sh
```
`lemmas.py` writes a cut-map cache `<circuit>.cm.pkl` next to the circuit; copy circuits to a work directory
first.

## Notes
- Cut map and signed substitution (paper, Section 5.2). `compositional_proof/cutmap.py` keeps one signed literal
  per word bit (the most-used simulation match; other matches stay ordinary interior nodes), and C2 of
  `check_compose.py` rejects an output bit produced twice, so no lemma has alternative literals.
  `lemmas.build_lemma` maps leaf nodes to free variables through a polarity-aware union-find (slots served by the
  same node share a variable), and in round 4 expands a_0 and e_0 through the candidate's logic down to the bits
  of W_0. The lemma builder is trusted: C7 reuses it, and C3/C4 compare unsigned node sets.
- The LRAT proofs are not included (up to 8.6 GB each). Manifests record their SHA-256; re-running the scripts
  regenerates proofs with the same CNF hashes (the proof bytes depend on the solver version).
- File names, environment variables and code comments were made neutral for publication; behaviour is
  unchanged (the generators reproduce the circuits byte for byte). Because comments were edited, the
  source-file hashes recorded inside the original proof manifests (`input_sha256.txt`, `tools` entries) refer
  to the original copies of the sources, not to the files shipped here; circuit, CNF and proof hashes are
  unaffected. The working-directory prefixes of the CNF and LRAT paths in the recorded solver commands were
  shortened to `work/`.
- The variable-CV circuit's own equivalence proof (an adapted compositional verifier for the 104,023-gate
  circuit, and a twin-cut certificate for the joint sigma pass) is not included; in the artifact it is checked
  by the GF(2) checkers and by simulation only.
