"""Run a protocol batch and retain OS resource use. No external dependencies.

Usage: python3 validation/uncertainty/run_batch.py evaluation 1 2 3
Omit scenario IDs to run every scenario. Working directory: repository root.
"""
import json
import platform
from pathlib import Path
import resource
import subprocess
import sys
import time

mode = sys.argv[1] if len(sys.argv) > 1 else "evaluation"
if mode not in ("evaluation", "smoke"):
    raise SystemExit("mode must be evaluation or smoke")
ids = [str(int(x)) for x in sys.argv[2:]]
out = Path("validation/uncertainty/results" if mode == "evaluation"
           else "validation/generated/uncertainty-smoke")
out.mkdir(parents=True, exist_ok=True)
command = ["Rscript", "validation/uncertainty/run_study.R", mode, *ids]
start = time.monotonic()
result = subprocess.run(command, check=False)
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
record = {
    "command": command, "mode": mode, "exit_code": result.returncode,
    "wall_seconds": time.monotonic() - start,
    "child_user_seconds": usage.ru_utime, "child_system_seconds": usage.ru_stime,
    "child_peak_rss": usage.ru_maxrss,
    "rss_units": "bytes" if platform.system() == "Darwin" else "KiB",
    "platform": platform.platform(), "python": platform.python_version(),
    "peak_scope": "maximum child-process RSS, not aggregate across concurrent processes",
}
name = "all" if not ids else "-".join(ids)
(out / f"resources-{name}.json").write_text(json.dumps(record, indent=2) + "\n")
raise SystemExit(result.returncode)
