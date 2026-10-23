"""Deterministic, unit-balanced sampling with rare-group coverage.

Fine labels are inclusion strata, never independent evaluation truth. This
module samples existing cells without replacement and does not infer patients.
"""
import hashlib
import numpy as np


def stable_seed(seed, source, pool):
    value = ('%s\0%s\0%s' % (seed, source, pool)).encode()
    return int.from_bytes(hashlib.sha256(value).digest()[:8], 'little')


def balanced_allocations(sizes, total, rng):
    sizes = np.asarray(sizes, dtype=np.int64)
    if sizes.ndim != 1 or (sizes < 0).any() or total < 0:
        raise ValueError('Invalid strata sizes or budget')
    total = min(int(total), int(sizes.sum()))
    allocated = np.zeros(len(sizes), dtype=np.int64)
    left = total
    while left:
        active = np.flatnonzero(allocated < sizes)
        if not len(active):
            raise ValueError('Allocation unexpectedly exhausted')
        active = rng.permutation(active)
        quota = max(1, left // len(active))
        for index in active:
            value = min(int(sizes[index] - allocated[index]), quota, left)
            allocated[index] += value
            left -= value
            if not left:
                break
    return allocated


def balanced_cell_selection(unit, fine, eligible, cap, seed, rare_reserve=2):
    unit = np.asarray(unit, dtype=object)
    fine = np.asarray(fine, dtype=object)
    eligible = np.asarray(eligible, dtype=bool)
    if unit.ndim != 1 or fine.shape != unit.shape or eligible.shape != unit.shape:
        raise ValueError('Unit/fine/eligibility axes differ')
    if cap < 0 or rare_reserve < 0:
        raise ValueError('Negative sampling budget')
    if any(str(x) in {'', 'nan', 'None'} for x in unit[eligible]):
        raise ValueError('Eligible cell lacks a declared sampling unit')
    rng = np.random.default_rng(seed)
    positions = np.flatnonzero(eligible)
    if not len(positions) or cap == 0:
        return np.array([], dtype=np.int64)
    names, unit_codes = np.unique(unit[positions].astype(str), return_inverse=True)
    sizes = np.bincount(unit_codes, minlength=len(names))
    budgets = balanced_allocations(sizes, cap, rng)
    selected = []
    for code, budget in enumerate(budgets):
        if not budget:
            continue
        rows = positions[unit_codes == code]
        _, fine_codes = np.unique(fine[rows].astype(str), return_inverse=True)
        counts = np.bincount(fine_codes)
        reserve_sizes = np.minimum(counts, rare_reserve)
        reserve = balanced_allocations(reserve_sizes, int(budget), rng)
        reserved = []
        for fine_code, amount in enumerate(reserve):
            if amount:
                candidates = rows[fine_codes == fine_code]
                reserved.extend(rng.choice(candidates, int(amount), replace=False).tolist())
        reserved = np.asarray(reserved, dtype=np.int64)
        remaining = int(budget) - len(reserved)
        if remaining:
            other = rows[~np.isin(rows, reserved)]
            reserved = np.concatenate([reserved, rng.choice(other, remaining, replace=False)])
        selected.extend(reserved.tolist())
    result = np.sort(np.asarray(selected, dtype=np.int64))
    if len(result) != min(int(cap), len(positions)) or len(np.unique(result)) != len(result):
        raise ValueError('Sampler changed budget or repeated an existing cell')
    if not eligible[result].all():
        raise ValueError('Sampler crossed the eligible population')
    return result
