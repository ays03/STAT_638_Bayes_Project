"""
Step 0: fetch the raw dataset.

The raw file is not tracked in git: it is 127 MB, above GitHub's 100 MB
per-file limit. This script downloads and unpacks it into data/ so the rest of
the pipeline can run from a fresh clone.

    python3 src/00_download.py
"""
import hashlib
import os
import sys
import urllib.request
import zipfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA = os.path.join(ROOT, "data")
URL = (
    "https://archive.ics.uci.edu/static/public/235/"
    "individual+household+electric+power+consumption.zip"
)
ZIP = os.path.join(DATA, "household_power_consumption.zip")
TXT = os.path.join(DATA, "household_power_consumption.txt")

# Expected shape of the raw file, as documented by UCI. Checked after download
# so a truncated or redirected fetch fails loudly here rather than silently
# producing a short series later.
EXPECTED_ROWS = 2_075_260  # 2,075,259 records plus the header
EXPECTED_BYTES = 132_960_755


def main():
    os.makedirs(DATA, exist_ok=True)

    if os.path.exists(TXT) and os.path.getsize(TXT) == EXPECTED_BYTES:
        print(f"already present and correct size: {TXT}")
        return

    if not os.path.exists(ZIP):
        print(f"downloading {URL}")
        try:
            with urllib.request.urlopen(URL, timeout=120) as r, open(ZIP, "wb") as f:
                total = 0
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
                    total += len(chunk)
                    print(f"\r  {total/1e6:7.1f} MB", end="", flush=True)
            print()
        except Exception as e:
            if os.path.exists(ZIP):
                os.remove(ZIP)
            sys.exit(f"download failed: {e}\nFetch it manually from "
                     "https://archive.ics.uci.edu/dataset/235/ into data/")

    print(f"unpacking {os.path.basename(ZIP)}")
    with zipfile.ZipFile(ZIP) as z:
        z.extractall(DATA)

    if not os.path.exists(TXT):
        sys.exit(f"expected {TXT} inside the archive but it is not there")

    size = os.path.getsize(TXT)
    with open(TXT, "rb") as f:
        n_rows = sum(1 for _ in f)
    print(f"  size {size:,} bytes (expected {EXPECTED_BYTES:,})")
    print(f"  rows {n_rows:,} (expected {EXPECTED_ROWS:,})")
    if size != EXPECTED_BYTES or n_rows != EXPECTED_ROWS:
        sys.exit("raw file does not match the documented dataset; stopping so "
                 "that downstream results are not silently wrong")
    print("\nok. Next: python3 src/01_aggregate.py")


if __name__ == "__main__":
    main()
