"""Same small-metadata scope, accept explicit English named calendar months.

No numeric day/month inference, two-digit century inference or fuzzy parsing.
The ISO-only failure remains historical; dates are not clinical acceptance.
"""
from datetime import date
from pathlib import Path
import hashlib
import importlib.util
import re

MONTHS = {name: i+1 for i, name in enumerate(
    ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'])}
PATH = Path(__file__).with_name('qualify_Hoch_melanoma_clinical_metadata_v1.py')
assert hashlib.sha256(PATH.read_bytes()).hexdigest() == '2013a1e450650b4630752c9d1cb3331dcde57b127d39b1b5b4832b97e65bd241'
SPEC = importlib.util.spec_from_file_location('frozen_Hoch_metadata_v1', PATH)
core = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(core)
ISO = core.iso_date


def explicit_date(value):
    if core.missing(value):return None
    text = str(value).strip()
    parsed = ISO(text)
    if parsed is not None:return parsed
    found = re.fullmatch(r'(\d{2})-([A-Za-z]{3})-(\d{4})', text)
    if not found or found[2].lower() not in MONTHS:return None
    try:return date(int(found[3]), MONTHS[found[2].lower()], int(found[1]))
    except ValueError:return None


def summarize(rows):
    previous = core.iso_date
    try:
        core.iso_date = explicit_date
        result = core.summarize(rows)
    finally:core.iso_date = previous
    for counts in result['source_date_format_counts'].values():
        if 'explicit_ISO_parsed' in counts:
            counts['explicit_ISO_or_named_calendar_month_parsed'] = counts.pop('explicit_ISO_parsed')
        if 'unparsed_nonISO_or_invalid' in counts:
            counts['unaccepted_calendar_representation_or_invalid'] = counts.pop('unparsed_nonISO_or_invalid')
    result.update(date_parser='Exact ISO or DD-EnglishMonthAbbreviation-YYYY; strict calendar validation',
                  no_numeric_day_month_or_century_guess=True,
                  conflicting_patient_contexts_are_unresolved_possible_multiple_episodes_not_bad_cases=True,
                  metadata_date_order_is_not_independent_clinical_endpoint_qualification=True,
                  date_wrapper_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    return result


if __name__ == '__main__':
    original = core.summarize
    # Keep the referenced original separate: summarize calls the frozen function.
    def selected(rows):
        core.summarize = original
        try:return summarize(rows)
        finally:core.summarize = selected
    core.summarize = selected
    core.main()
