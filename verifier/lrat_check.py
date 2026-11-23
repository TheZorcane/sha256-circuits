"""Strict, independent, RUP-only LRAT proof checker (text LRAT).

This checker is written independently of the SAT solver and of any other
proof checker; it is the trust anchor for every "equivalent" verdict.  The
harness (equiv.py) additionally runs a second, external checker.

Rules (deliberately narrow; CaDiCaL --lrat emits only RUP steps):
  * addition line: "<id> <lits...> 0 <hints...> 0", id strictly greater than
    every clause id seen so far (input clauses are 1..m);
  * all hints must be positive ids of live clauses (negative hints = RAT steps
    are REJECTED);
  * with the lemma's literals assumed false, each hint clause must, in order,
    be unit (its single unassigned literal becomes true) or falsified
    (conflict: lemma accepted).  A satisfied hint, a hint with >= 2
    unassigned literals, or running out of hints without conflict => reject;
  * deletion line: "<id> d <ids...> 0"; ids must be live;
  * the proof is accepted iff the empty clause is added and no line failed.
Soundness: each accepted lemma is implied by the live clauses (reverse unit
propagation), so the empty clause is implied by the CNF => CNF unsatisfiable.
"""

import sys


class ProofError(Exception):
    pass


def parse_cnf(path):
    clauses, nv, nc = [], None, None
    cur = []
    with open(path) as f:
        for ln in f:
            if ln.startswith("c"):
                continue
            if ln.startswith("p"):
                t = ln.split()
                if t[:2] != ["p", "cnf"]:
                    raise ProofError("bad CNF header")
                nv, nc = int(t[2]), int(t[3])
                continue
            for tok in ln.split():
                v = int(tok)
                if v == 0:
                    clauses.append(cur)
                    cur = []
                else:
                    cur.append(v)
    if cur or nc is None or len(clauses) != nc:
        raise ProofError("CNF clause count mismatch")
    return clauses


def check(cnf_path, proof_path):
    """Return (accepted: bool, message: str)."""
    try:
        db = {i + 1: c for i, c in enumerate(parse_cnf(cnf_path))}
        last_id = len(db)
        n_added = n_deleted = 0
        with open(proof_path) as f:
            for lineno, ln in enumerate(f, 1):
                t = ln.split()
                if not t or t[0] == "c":
                    continue
                try:
                    cid = int(t[0])
                except ValueError:
                    raise ProofError(f"line {lineno}: bad id")
                if len(t) > 1 and t[1] == "d":
                    ids = [int(x) for x in t[2:]]
                    if not ids or ids[-1] != 0:
                        raise ProofError(f"line {lineno}: deletion not 0-terminated")
                    for d in ids[:-1]:
                        if d not in db:
                            raise ProofError(f"line {lineno}: deleting non-live clause {d}")
                        del db[d]
                        n_deleted += 1
                    continue
                nums = [int(x) for x in t[1:]]
                if nums.count(0) != 2 or nums[-1] != 0:
                    raise ProofError(f"line {lineno}: malformed addition")
                z = nums.index(0)
                lemma, hints = nums[:z], nums[z + 1:-1]
                if cid <= last_id:
                    raise ProofError(f"line {lineno}: id {cid} not increasing")
                true = {-l for l in lemma}           # literals currently true
                conflict = False
                for h in hints:
                    if h <= 0:
                        raise ProofError(f"line {lineno}: non-positive hint {h} (RAT not supported)")
                    c = db.get(h)
                    if c is None:
                        raise ProofError(f"line {lineno}: hint {h} not a live clause")
                    unassigned = None
                    for l in c:
                        if -l in true:
                            continue                 # literal false
                        if l in true:
                            raise ProofError(f"line {lineno}: hint {h} already satisfied")
                        if unassigned is not None and unassigned != l:
                            raise ProofError(f"line {lineno}: hint {h} not unit")
                        unassigned = l
                    if unassigned is None:
                        conflict = True
                        break
                    true.add(unassigned)
                if not conflict:
                    raise ProofError(f"line {lineno}: no conflict for lemma {cid}")
                db[cid] = lemma
                last_id = cid
                n_added += 1
                if not lemma:
                    return True, f"VERIFIED: empty clause {cid} derived ({n_added} lemmas, {n_deleted} deletions checked)"
        return False, "NOT VERIFIED: empty clause never derived"
    except ProofError as e:
        return False, f"NOT VERIFIED: {e}"
    except (ValueError, OSError) as e:
        return False, f"NOT VERIFIED: {type(e).__name__}: {e}"


if __name__ == "__main__":
    ok, msg = check(sys.argv[1], sys.argv[2])
    print(("s VERIFIED " if ok else "s NOT VERIFIED ") + msg)
    sys.exit(0 if ok else 1)
