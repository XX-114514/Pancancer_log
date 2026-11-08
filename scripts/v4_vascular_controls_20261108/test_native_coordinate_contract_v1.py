"""Regression checks for the actual author unit mismatch; no human data."""
from types import SimpleNamespace
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from native_coordinate_contract_v1 import native_micron_coordinates

xy = np.array([[100., 200.], [145., 200.], [200., 260.]])


def obj(values):
    return SimpleNamespace(
        obs=pd.DataFrame(xy, columns=['x_centroid', 'y_centroid']),
        obsm={'xenium_um_spatial': values, 'spatial': xy * 7 + 1234},
    )


native, receipt = native_micron_coordinates(obj(xy / .2125))
assert np.array_equal(native, xy)
assert receipt['archived_actual_units'] == 'native_pixels_despite_um_key_name'
assert receipt['maximum_absolute_crosscheck_error_um'] < 1e-10
assert cKDTree(native).query_pairs(50) == {(0, 1)}
assert cKDTree(xy / .2125).query_pairs(50) == set()
assert native_micron_coordinates(obj(xy))[1]['archived_to_micrometres_multiplier'] == 1
for bad in [xy * 2, xy[:, ::-1], xy[::-1], xy + [10, 20], xy[:2], xy * np.nan]:
    try:
        native_micron_coordinates(obj(bad))
    except ValueError:
        continue
    raise AssertionError('Unknown or invalid coordinates silently accepted')
missing = obj(xy)
missing.obsm.pop('xenium_um_spatial')
try:
    native_micron_coordinates(missing)
except ValueError:
    pass
else:
    raise AssertionError('Missing crosscheck silently accepted')
for invalid_scale in [0, -1, np.nan]:
    try:
        native_micron_coordinates(obj(xy), author_pixel_size_um=invalid_scale)
    except ValueError:
        continue
    raise AssertionError('Invalid documented scale silently accepted')
print('PASS_SYNTHETIC_ONLY: author pixels, genuine microns, 50um edges, image separation; unknown scale/axis/order/translation/NaN/shape/missing metadata rejected')
