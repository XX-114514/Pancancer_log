#!/usr/bin/env python3
"""Emit download-manifest rows (http_code, bytes, relative_dir, url) for one GEO series.

Reads the series filelist.txt (per-GSM supplementary files with sizes). Files are kept when they
match --include and not --exclude; the RAW.tar archive is never listed (individual files instead).
If the series has no filelist (series-level files only), --series-files lists those with size 0.
Stdlib only.  Usage:
  python3 -I make_geo_manifest.py GSE299786 S20_GSE299786_CosMx --exclude '_tx_file' >> manifest.tsv
"""
import argparse
import re
import sys
import urllib.error
import urllib.request

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"


def stub(acc):
    return acc[:-3] + "nnn"


def fetch(url):
    with urllib.request.urlopen(url, timeout=60) as fh:
        return fh.read().decode("utf-8", "replace")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gse")
    ap.add_argument("dest_prefix", help="top-level directory under RAW_ROOT")
    ap.add_argument("--include", default=".")
    ap.add_argument("--exclude", default=None)
    ap.add_argument("--series-files", action="store_true", help="list series-level suppl files instead")
    a = ap.parse_args()
    inc = re.compile(a.include)
    exc = re.compile(a.exclude) if a.exclude else None
    keep = lambda n: inc.search(n) and not (exc and exc.search(n))
    base = f"{GEO}/series/{stub(a.gse)}/{a.gse}/suppl/"
    rows = []
    if a.series_files:
        for name in sorted(set(re.findall(r'href="([^"/?]+\.[A-Za-z0-9.]+)"', fetch(base)))):
            if keep(name) and name != "filelist.txt" and not name.endswith("_RAW.tar"):
                rows.append(("200", "0", a.dest_prefix, base + name))
    else:
        try:
            listing = fetch(base + "filelist.txt")
        except urllib.error.HTTPError as e:
            sys.exit(f"{a.gse}: no filelist.txt ({e.code}); rerun with --series-files")
        for line in listing.splitlines()[1:]:
            parts = line.rstrip("\r").split("\t")
            if len(parts) < 4 or parts[0] != "File" or not keep(parts[1]):
                continue
            name, size = parts[1], parts[3]
            gsm = name.split("_", 1)[0]
            if not gsm.startswith("GSM"):
                continue
            url = f"{GEO}/samples/{stub(gsm)}/{gsm}/suppl/{name}"
            rows.append(("200", size, f"{a.dest_prefix}/{gsm}", url))
    if not rows:
        sys.exit(f"{a.gse}: no files matched")
    for r in rows:
        print("\t".join(r))
    print(f"{a.gse}: {len(rows)} files", file=sys.stderr)


if __name__ == "__main__":
    main()
