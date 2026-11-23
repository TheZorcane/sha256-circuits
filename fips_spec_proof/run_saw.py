"""Run one SAW script with a wall-clock budget; on timeout kill SAW and its solver children (own process
group only).  usage: run_saw.py TIMEOUT_S LOGFILE SCRIPT.saw
Exit code: SAW's exit code, or 124 on timeout or when the group RSS exceeds SAW_MEM_MB (default 16000).  CRYPTOLPATH must be set by the caller."""
import os
import signal
import subprocess
import sys
import time

SAW = os.environ.get("SAW", "saw")


MEM_KB = int(os.environ.get("SAW_MEM_MB", "16000")) * 1024


def group_rss_kb(pgid):
    out = subprocess.run(["ps", "-A", "-o", "pgid=,rss="], capture_output=True, text=True).stdout
    return sum(int(r) for g, r in (ln.split() for ln in out.splitlines() if ln.strip()) if int(g) == pgid)


def main():
    budget, log, script = float(sys.argv[1]), sys.argv[2], sys.argv[3]
    t0 = time.time()
    with open(log, "w") as f:
        f.write(f"# run_saw: {script} budget={budget:.0f}s start={time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.flush()
        p = subprocess.Popen([SAW, os.path.basename(script)], cwd=os.path.dirname(os.path.abspath(script)),
                             stdout=f, stderr=subprocess.STDOUT, start_new_session=True)
        rc, why = None, ""
        while rc is None:
            try:
                rc = p.wait(timeout=10)
            except subprocess.TimeoutExpired:
                rss = group_rss_kb(p.pid)
                if time.time() - t0 > budget or rss > MEM_KB:
                    why = "TIMEOUT" if time.time() - t0 > budget else f"MEMORY LIMIT ({rss // 1024} MB)"
                    os.killpg(p.pid, signal.SIGKILL)
                    p.wait()
                    rc = 124
        f.write(f"\n# run_saw: exit={rc} wall={time.time() - t0:.1f}s{' (' + why + ')' if why else ''}\n")
    print(f"{os.path.basename(script)}: exit={rc} wall={time.time() - t0:.1f}s")
    sys.exit(rc)


if __name__ == "__main__":
    main()
