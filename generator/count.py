"""count.py F.xag...: verifier/cost.report_xag of the canonical form (counts, depth, sha256)."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "verifier"))
import cost
from xag import read_xag
for p in sys.argv[1:]:
    x = read_xag(open(p).read()); c = cost.canonicalize(x); r = cost.report_xag(c)
    r["canonical_input"] = cost.serialize(c) == open(p, "rb").read()
    print(os.path.basename(p), json.dumps(r))
