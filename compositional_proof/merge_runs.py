"""merge_runs.py DST SRC1 SRC2 ...: merge lemma runs of the SAME candidate and verifier that were started with
different --only lists: copy every per-lemma directory into DST and concatenate the lemmas.json records
(duplicates are an error).  check_compose.py then re-checks every record against its manifest."""
import sys, os, json, shutil
dst = sys.argv[1]; os.makedirs(dst, exist_ok=True)
recs = {}
for src in sys.argv[2:]:
    for r in json.load(open(os.path.join(src, "lemmas.json"))):
        assert r["lemma"] not in recs, r["lemma"]
        recs[r["lemma"]] = r
        d = r["lemma"].replace(":", "_")
        if not os.path.exists(os.path.join(dst, d)):
            shutil.copytree(os.path.join(src, d), os.path.join(dst, d))
json.dump([recs[k] for k in sorted(recs)], open(os.path.join(dst, "lemmas.json"), "w"))
print(len(recs), "lemmas merged into", dst)
