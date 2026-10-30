"""Readable tables from the complete CSV summaries; no scenario selection."""
import csv
import json
from pathlib import Path

root = Path(__file__).resolve().parent
p = json.loads((root / "protocol-v1.json").read_text())
rows = list(csv.DictReader((root / "results" / "summary.csv").open()))
lookup = {(r["scenario"],r["method"]):r for r in rows}
text = ["# Complete uncertainty study tables", "",
        "Percentages with Wilson 95% Monte Carlo intervals. These describe the specified simulations, not universal performance.", "",
        "## Primary-target interval coverage", "",
        "| Scenario | Locus | Block | Individual | Crossed |", "|---|---:|---:|---:|---:|"]
for s in p["scenarios"]:
    cells = []
    for m in p["methods"]:
        r=lookup[(s["name"],m)]
        cells.append(f"{100*float(r['coverage']):.2f} [{100*float(r['coverage_lower']):.2f}, {100*float(r['coverage_upper']):.2f}]")
    text.append("| "+s["name"]+" | "+" | ".join(cells)+" |")
text += ["", "Fixed panels use the realized-panel target; random panels use the genomic-distribution target. See summary.csv for secondary panel coverage.", "",
         "## False genetic-direction support under the null", "", "| Scenario | Locus | Block | Individual | Crossed |", "|---|---:|---:|---:|---:|"]
for s in p["scenarios"]:
    if not s["name"].startswith("null_"): continue
    cells=[]
    for m in p["methods"]:
        r=lookup[(s["name"],m)]
        cells.append(f"{100*float(r['supported_rate']):.2f} [{100*float(r['supported_lower']):.2f}, {100*float(r['supported_upper']):.2f}]")
    text.append("| "+s["name"]+" | "+" | ".join(cells)+" |")
text += ["", "Each null has 1,000 datasets. The null is zero genetic contrast; demographic migration is not simulated.", "",
         "## Alternative direction counts", "", "Cells are correct / wrong / unresolved among 625 datasets.", "",
         "| Scenario | Locus | Block | Individual | Crossed |", "|---|---:|---:|---:|---:|"]
for s in p["scenarios"]:
    if s["name"].startswith("null_"): continue
    cells=[]
    for m in p["methods"]:
        r=lookup[(s["name"],m)]
        cells.append(" / ".join(r[k] for k in ("correct","wrong","unresolved")))
    text.append("| "+s["name"]+" | "+" | ".join(cells)+" |")
text += ["", "Perfect sign recovery in these alternatives coexists with poor interval coverage. It does not validate migration inference or general performance.", ""]
(root / "results" / "tables.md").write_text("\n".join(text))
