"""Archive all planned outcomes with reproducible ZIP metadata and checksums.

No row, scenario, failure, or negative outcome is filtered here.
"""
import csv
import hashlib
import json
from pathlib import Path
import zipfile

root = Path(__file__).resolve().parent
out = root / "results"
protocol = json.loads((root / "protocol-v1.json").read_text())
completion = json.loads((out / "completion.json").read_text())
assert completion["complete"] and completion["datasets"] == 15125
assert completion["method_results"] == 60500 and completion["execution_failures"] == 0
archive = out / "replicates.zip"
records = []
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for scenario in protocol["scenarios"]:
        path = out / (scenario["name"] + ".csv")
        data = path.read_bytes()
        rows = list(csv.DictReader(data.decode().splitlines()))
        n = protocol["null_replicates"] if scenario["name"].startswith("null_") else protocol["alternative_replicates"]
        assert len(rows) == n * len(protocol["methods"])
        assert all(r["status"] == "ok" for r in rows)
        info = zipfile.ZipInfo(path.name, date_time=(1980, 1, 1, 0, 0, 0))
        info.compress_type = zipfile.ZIP_DEFLATED
        info.create_system = 3
        info.external_attr = 0o100644 << 16
        z.writestr(info, data, compresslevel=9)
        records.append({"file": path.name, "rows": len(rows), "bytes": len(data),
                        "sha256": hashlib.sha256(data).hexdigest()})
with zipfile.ZipFile(archive) as z:
    assert len(z.namelist()) == len(records)
    for record in records:
        assert hashlib.sha256(z.read(record["file"])).hexdigest() == record["sha256"]
manifest = {"archive": archive.name, "bytes": archive.stat().st_size,
            "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
            "datasets": 15125, "method_results": 60500, "files": records}
(out / "replicate-archive-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(f"Verified {len(records)} CSVs, 60,500 method rows; archive {archive.stat().st_size:,} bytes.")
