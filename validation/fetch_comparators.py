"""Fetch exact upstream source files for an optional, isolated numeric audit.

Run from the repository root. No upstream code is vendored or shipped in the
R package. Original sources retain their upstream licenses. Network downloads
are accepted only when their Git blob IDs match the committed manifest.
"""

import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parent


def blob_id(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def main():
    pins = json.loads((ROOT / "comparator-pins.json").read_text())
    records = []
    for pin in pins:
        target = ROOT / "generated" / "comparators" / pin["software"] / pin["path"]
        url = (f"https://raw.githubusercontent.com/{pin['repository']}/"
               f"{pin['commit']}/{pin['path']}")
        if target.exists():
            data = target.read_bytes()
            if blob_id(data) != pin["git_blob_sha"]:
                raise RuntimeError(f"Cached source checksum mismatch: {target}")
        else:
            request = urllib.request.Request(url, headers={"User-Agent": "snpflowR-comparator-audit"})
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read()
            if blob_id(data) != pin["git_blob_sha"]:
                raise RuntimeError(f"Downloaded source checksum mismatch: {url}")
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        records.append({**pin, "url": url, "sha256": hashlib.sha256(data).hexdigest()})
    destination = ROOT / "results" / "comparator-source-manifest.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(records, indent=2) + "\n")
    print(f"Verified {len(records)} upstream source blobs.")


if __name__ == "__main__":
    main()
