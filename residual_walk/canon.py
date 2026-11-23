"""canon.py IN.xag OUT.xag : verifier/cost.py canonicalize + report (no other rewriting)."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
import cost, xag
x = xag.read_xag(open(sys.argv[1]).read())
y = cost.canonicalize(x)
open(sys.argv[2], "w").write(xag.write_xag(y))
print(json.dumps(cost.report_xag(y)))
