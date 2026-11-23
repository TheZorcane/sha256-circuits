"""Combinational equivalence harness: miter -> Tseitin CNF -> CaDiCaL with
proof -> independent proof check -> JSON manifest.

Interface map: both circuits are first brought to a common interface with
sim.normalize(x, in_map, out_map, n) (or are already in it).  The miter shares
the n inputs and asserts OR_i (A_i XOR B_i): UNSAT <=> A == B on all inputs.

Status values
  equivalent        UNSAT and the LRAT proof was accepted by BOTH verifier/lrat_check.py
                    (strict RUP) and the vendored lrat-check
  unsat_unchecked   UNSAT but no proof check succeeded
  not_equivalent    SAT and the counterexample was confirmed by simulation
  sat_unconfirmed   SAT but simulation did NOT confirm (harness bug: investigate)
  unverified        timeout / solver error
Trust boundary for "equivalent": this file's miter + Tseitin encoding, the
XAG translation of each circuit (xag.py), and lrat_check.py (the strict,
independent RUP checker).  The vendored lrat-check is a second opinion only.
The solver is NOT trusted.
"""

import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import time

import cost
import lrat_check
import sim
from xag import AND, XAG

HERE = os.path.dirname(os.path.abspath(__file__))
VENDOR = os.path.join(HERE, "vendor")
CADICAL = os.environ.get("CADICAL", "cadical")
LRAT_CHECK = os.path.join(VENDOR, "lrat-check")
DRAT_TRIM = os.path.join(VENDOR, "drat-trim")


# --------------------------------------------------------------------------- #
# miter + Tseitin
# --------------------------------------------------------------------------- #
def miter_cnf(a: XAG, b: XAG):
    """Return (dimacs_text, input_vars) where input_vars[i] is the CNF variable
    of shared input i.  CNF var 1 is constant FALSE."""
    if a.n_inputs != b.n_inputs or len(a.outputs) != len(b.outputs):
        raise ValueError("interface mismatch: normalise both circuits first")
    n = a.n_inputs
    clauses = [[-1]]                    # var 1 = FALSE
    inputs = list(range(2, 2 + n))
    nxt = [2 + n]

    def encode(x):
        var = [1] + inputs[:]           # XAG var -> CNF var

        def lit(l):
            v = var[l >> 1]
            return -v if l & 1 else v

        for op, p, q in x.nodes:
            o = nxt[0]; nxt[0] += 1
            la, lb = lit(p), lit(q)
            if op == AND:
                clauses.extend([[-o, la], [-o, lb], [o, -la, -lb]])
            else:
                clauses.extend([[-o, la, lb], [-o, -la, -lb],
                                [o, -la, lb], [o, la, -lb]])
            var.append(o)
        return [lit(l) for l in x.outputs]

    oa, ob = encode(a), encode(b)
    diffs = []
    for la, lb in zip(oa, ob):
        d = nxt[0]; nxt[0] += 1
        clauses.extend([[-d, la, lb], [-d, -la, -lb], [d, -la, lb], [d, la, -lb]])
        diffs.append(d)
    clauses.append(diffs)
    nv = nxt[0] - 1
    out = [f"p cnf {nv} {len(clauses)}"]
    out += [" ".join(map(str, c)) + " 0" for c in clauses]
    return "\n".join(out) + "\n", inputs


# --------------------------------------------------------------------------- #
# process runner with per-child rusage (peak RSS)
# --------------------------------------------------------------------------- #
REPO = os.path.dirname(HERE)


def _rel(p):
    """Repo-relative form of absolute paths inside the repo (for manifests)."""
    if os.path.isabs(p) and p.startswith(REPO + os.sep):
        return os.path.relpath(p, REPO)
    return p


def run(cmd, timeout, stdout_path=None):
    t0 = time.time()
    out = open(stdout_path, "w") if stdout_path else subprocess.DEVNULL
    p = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT)
    timed_out = False
    while True:
        pid, status, ru = os.wait4(p.pid, os.WNOHANG)
        if pid:
            break
        if time.time() - t0 > timeout:
            p.kill()
            pid, status, ru = os.wait4(p.pid, 0)
            timed_out = True
            break
        time.sleep(0.02)
    if stdout_path:
        out.close()
    rss = ru.ru_maxrss if platform.system() == "Darwin" else ru.ru_maxrss * 1024
    code = os.waitstatus_to_exitcode(status)
    shown = [os.path.basename(cmd[0])] + [_rel(c) for c in cmd[1:]]
    return {"cmd": shown, "exit_code": code,
            "timed_out": timed_out, "wall_s": round(time.time() - t0, 3),
            "peak_rss_bytes": rss}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def tool_versions():
    v = {}
    try:
        v["cadical"] = subprocess.run([CADICAL, "--version"], capture_output=True,
                                      text=True).stdout.strip()
    except OSError:
        v["cadical"] = None
    commit = open(os.path.join(VENDOR, "DRAT_TRIM_COMMIT.txt")).read().strip()
    for name, path in (("lrat-check", LRAT_CHECK), ("drat-trim", DRAT_TRIM)):
        v[name] = {"used": name == "lrat-check", "source_commit": commit,
                   "binary_sha256": sha256_file(path) if os.path.exists(path) else None,
                   "source_sha256": sha256_file(os.path.join(VENDOR, name + ".c"))}
    v["python"] = platform.python_version()
    return v


def check_lrat(cnf, proof, timeout, log):
    """Run the vendored lrat-check (second opinion).  Its exit code is not
    used.  Accepted only if some
    output line is exactly 'c VERIFIED', no line mentions failure or
    'NOT VERIFIED', and the run did not time out.  The proof must be text LRAT
    (cadical is run with --binary=false)."""
    chk = run([LRAT_CHECK, cnf, proof], timeout, stdout_path=log)
    text = open(log).read()
    lines = [l.strip() for l in text.splitlines()]
    ok = (not chk["timed_out"]
          and "c VERIFIED" in lines
          and not any("FAIL" in l.upper() or "NOT VERIFIED" in l.upper() for l in lines))
    chk["verdict_line"] = next((l for l in text.splitlines()
                                if "VERIFIED" in l or "FAIL" in l.upper()), None)
    return dict(chk, checker="lrat-check", format="LRAT (text)", accepted=ok)


def check_proof(cnf, proof, timeout, log):
    """Proof accepted only if BOTH the independent strict RUP checker
    (lrat_check.py) and the vendored lrat-check accept.  The strict checker is
    the trust anchor; the vendored checker is a second opinion, never
    sufficient by itself."""
    t0 = time.time()
    py_ok, py_msg = lrat_check.check(cnf, proof)
    py = {"checker": "verifier/lrat_check.py (strict RUP-only LRAT)", "accepted": py_ok,
          "message": py_msg, "wall_s": round(time.time() - t0, 3)}
    c = check_lrat(cnf, proof, timeout, log)
    return {"accepted": py_ok and c["accepted"], "rule": "both checkers must accept",
            "python_checker": py, "vendored_lrat_check": c}


# --------------------------------------------------------------------------- #
# main entry
# --------------------------------------------------------------------------- #
def check(a: XAG, b: XAG, workdir, name="miter", timeout=600, check_timeout=None,
          interface=None, extra=None):
    """Prove/refute a == b (same interface).  Writes CNF/proof/logs to workdir,
    returns the manifest dict (also written to workdir/<name>.manifest.json)."""
    os.makedirs(workdir, exist_ok=True)
    check_timeout = check_timeout or max(4 * timeout, 60)
    t0 = time.time()
    cnf_text, inputs = miter_cnf(a, b)
    cnf = os.path.join(workdir, name + ".cnf")
    with open(cnf, "w") as f:
        f.write(cnf_text)
    proof = os.path.join(workdir, name + ".lrat")
    log = os.path.join(workdir, name + ".cadical.log")
    m = {
        "name": name,
        "circuit_a_sha256": cost.xag_hash(a), "circuit_b_sha256": cost.xag_hash(b),
        "circuit_a_counts": a.counts(), "circuit_b_counts": b.counts(),
        "interface": interface or "identity (both circuits already in the same interface)",
        "cnf_sha256": hashlib.sha256(cnf_text.encode()).hexdigest(),
        "cnf_vars_clauses": cnf_text.split("\n", 1)[0],
        "tools": tool_versions(),
        "timeout_s": timeout,
        "host": f"{platform.system()} {platform.machine()}",
    }
    if extra:
        m.update(extra)
    solve = run([CADICAL, "--lrat=true", "--binary=false", "-q", cnf, proof],
                timeout, stdout_path=log)
    m["solver_run"] = solve
    if solve["timed_out"] or solve["exit_code"] not in (10, 20):
        m["result"] = "unknown"
        m["status"] = "unverified"
    elif solve["exit_code"] == 10:
        m["result"] = "sat"
        model = set()
        for ln in open(log):
            if ln.startswith("v "):
                model.update(int(t) for t in ln.split()[1:])
        bits = [1 if v in model else 0 for v in inputs]
        oa, ob = sim.eval_bits(a, bits), sim.eval_bits(b, bits)
        diff = [i for i in range(len(oa)) if oa[i] != ob[i]]
        m["counterexample_inputs"] = "".join(map(str, bits))
        m["counterexample_differing_outputs"] = diff
        m["status"] = "not_equivalent" if diff else "sat_unconfirmed"
    else:
        m["result"] = "unsat"
        m["proof_sha256"] = sha256_file(proof)
        m["proof_bytes"] = os.path.getsize(proof)
        m["proof_check"] = check_proof(cnf, proof, check_timeout,
                                       os.path.join(workdir, name + ".lrat-check.log"))
        m["status"] = "equivalent" if m["proof_check"]["accepted"] else "unsat_unchecked"
    m["total_wall_s"] = round(time.time() - t0, 3)
    with open(os.path.join(workdir, name + ".manifest.json"), "w") as f:
        json.dump(m, f, indent=2)
    return m


if __name__ == "__main__":
    import argparse
    from xag import read_xag
    ap = argparse.ArgumentParser(description="prove two .xag files equivalent "
                                 "(identical interfaces assumed)")
    ap.add_argument("a"); ap.add_argument("b"); ap.add_argument("workdir")
    ap.add_argument("--timeout", type=float, default=600)
    args = ap.parse_args()
    A = read_xag(open(args.a).read()); B = read_xag(open(args.b).read())
    r = check(A, B, args.workdir, timeout=args.timeout)
    print(json.dumps({k: r[k] for k in ("status", "result", "total_wall_s")}))
    sys.exit(0 if r["status"] in ("equivalent", "not_equivalent") else 1)
