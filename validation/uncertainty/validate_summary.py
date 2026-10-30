"""Independent-language checks of archived counts, coverage and Wilson limits."""
import csv
import io
import json
import math
from pathlib import Path
import statistics
import zipfile

root = Path(__file__).resolve().parent
out = root / "results"
p = json.loads((root / "protocol-v1.json").read_text())
summary = {(r["scenario"], r["method"]): r for r in csv.DictReader((out / "summary.csv").open())}
checked = 0
largest = 0.0
z = statistics.NormalDist().inv_cdf(.975)

def wilson(k, n):
    center = (k/n + z*z/(2*n))/(1+z*z/n)
    half = z*math.sqrt((k/n)*(1-k/n)/n + z*z/(4*n*n))/(1+z*z/n)
    return center-half, center+half

with zipfile.ZipFile(out / "replicates.zip") as archive:
    for s in p["scenarios"]:
        rows = list(csv.DictReader(io.StringIO(archive.read(s["name"]+".csv").decode())))
        seen = set()
        for r in rows:
            key = (int(r["replicate"]), r["method"])
            assert key not in seen
            seen.add(key)
            assert int(r["simulation_seed"]) == p["evaluation_seed_base"]+p["scenario_seed_stride"]*s["id"]+key[0]
            assert int(r["bootstrap_seed"]) == p["bootstrap_seed_base"]+p["scenario_seed_stride"]*s["id"]+key[0]
            lo, hi, truth, panel = (float(r[x]) for x in ("lower","upper","truth","panel_truth"))
            assert lo <= hi and all(math.isfinite(x) for x in (lo,hi,truth,panel))
            assert (r["covered"]=="TRUE") == (lo <= truth <= hi)
            assert (r["panel_covered"]=="TRUE") == (lo <= panel <= hi)
            expected = 0
            if r["eligible"]=="TRUE":
                if lo > 0 and float(r["estimate"]) > 0: expected = 1
                if hi < 0 and float(r["estimate"]) < 0: expected = -1
            assert int(r["direction"]) == expected
            checked += 1
        for method in p["methods"]:
            d = [r for r in rows if r["method"]==method]
            q = summary[(s["name"],method)]
            n = len(d)
            coverage = sum(r["covered"]=="TRUE" for r in d)
            support = sum(int(r["direction"])!=0 for r in d)
            assert n == int(q["n"]) and coverage == int(q["covered"]) and support == int(q["supported"])
            for field, k in (("coverage",coverage),("supported",support)):
                lower, upper = wilson(k,n)
                error = max(abs(lower-float(q[field+"_lower"])),abs(upper-float(q[field+"_upper"])))
                largest = max(largest,error)
                assert error < 1e-12
assert checked == 60500
(out / "summary-verification.json").write_text(json.dumps({
    "method_rows_checked": checked, "unique_scenarios": len(p["scenarios"]),
    "seed_and_interval_checks": "passed", "summary_counts": "passed",
    "maximum_wilson_difference": largest,
    "scope": "separate Python arithmetic check, not external scientific review"
}, indent=2)+"\n")
print(f"Verified {checked:,} archived method outcomes and all Wilson summary limits.")
