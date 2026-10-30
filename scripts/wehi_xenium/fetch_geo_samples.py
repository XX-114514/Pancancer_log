#!/usr/bin/env python3
"""Write one row per GSM (title, patient, run, source, characteristics) for a GEO series.

Stdlib only. GEO text responses use CRLF; values are stripped before use.
Patient/run are parsed from titles shaped like "Patient 9, run 2"; other titles leave them empty.
Usage: python3 -I fetch_geo_samples.py GSE300147 OUT_TSV
"""
import csv
import re
import sys
import urllib.request

URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ=gsm&form=text&view=brief"
TITLE = re.compile(r"^(?P<patient>.+?),\s*run\s*(?P<run>\w+)$", re.I)


def parse(text):
    rows, cur = [], None
    for line in text.splitlines():
        line = line.rstrip("\r")
        key, _, value = line.partition(" = ")
        value = value.strip()
        if key == "^SAMPLE":
            cur = {"gsm": value, "title": "", "patient": "", "run": "", "source": ""}
            rows.append(cur)
        elif cur is None:
            continue
        elif key == "!Sample_title":
            cur["title"] = value
            m = TITLE.match(value)
            if m:
                cur["patient"], cur["run"] = m.group("patient"), m.group("run")
        elif key == "!Sample_source_name_ch1":
            cur["source"] = value
        elif key == "!Sample_characteristics_ch1" and ": " in value:
            k, v = value.split(": ", 1)
            cur["char_" + re.sub(r"\W+", "_", k.strip().lower())] = v.strip()
    return rows


def main(gse, out_tsv):
    with urllib.request.urlopen(URL.format(gse), timeout=60) as fh:
        rows = parse(fh.read().decode("utf-8", "replace"))
    if not rows:
        sys.exit(f"no GSM records returned for {gse}")
    fields = ["gsm", "title", "patient", "run", "source"]
    fields += sorted({k for r in rows for k in r} - set(fields))
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"{gse}: {len(rows)} samples -> {out_tsv}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
