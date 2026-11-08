"""Coordinate contract for the P08 breast Xenium author objects.

The author's image preprocessing divides native microns by 0.2125 before
saving the misleadingly named ``xenium_um_spatial`` key. Never infer units
from that name or use the aligned H&E ``spatial`` pixels for distances.
"""
import numpy as np


def native_micron_coordinates(adata, author_pixel_size_um=0.2125, atol_um=1e-4):
    """Return unchanged native obs microns, after an explicit unit check.

    Accept either genuinely micron-valued legacy coordinates or the known
    author-native pixel representation. Unknown transformations fail closed;
    no fitted scale, translation, rotation, or row realignment is applied.
    """
    if not np.isfinite(author_pixel_size_um) or author_pixel_size_um <= 0:
        raise ValueError('A positive documented author pixel size is required')
    if not np.isfinite(atol_um) or atol_um <= 0:
        raise ValueError('A positive absolute micron tolerance is required')
    coords = adata.obs[['x_centroid', 'y_centroid']].to_numpy(dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 2 or not len(coords):
        raise ValueError('Native obs coordinates must be a nonempty n-by-2 array')
    if not np.isfinite(coords).all() or not np.all(np.ptp(coords, axis=0) > 0):
        raise ValueError('Native obs coordinates are nonfinite or degenerate')
    if 'xenium_um_spatial' not in adata.obsm:
        raise ValueError('Required author coordinate cross-check is unavailable')
    archived = np.asarray(adata.obsm['xenium_um_spatial'], dtype=float)
    if archived.shape != coords.shape or not np.isfinite(archived).all():
        raise ValueError('Author coordinate shape or finite-value check failed')
    direct_error = float(np.max(np.abs(coords - archived)))
    restored_error = float(np.max(np.abs(coords - archived * author_pixel_size_um)))
    if direct_error <= atol_um:
        archived_units = 'native_micrometres'
        archived_to_um = 1.0
        accepted_error = direct_error
    elif restored_error <= atol_um:
        archived_units = 'native_pixels_despite_um_key_name'
        archived_to_um = float(author_pixel_size_um)
        accepted_error = restored_error
    else:
        raise ValueError('Unrecognized author coordinate transform; reject input')
    return coords, {
        'status': 'PASS_native_obs_microns_with_explicit_author_unit_crosscheck',
        'distance_coordinate_source': 'obs[x_centroid,y_centroid]',
        'distance_coordinate_units': 'micrometres',
        'distance_coordinates_modified': False,
        'archived_key': 'xenium_um_spatial',
        'archived_actual_units': archived_units,
        'archived_to_micrometres_multiplier': archived_to_um,
        'author_documented_micrometres_per_pixel': float(author_pixel_size_um),
        'absolute_tolerance_um': float(atol_um),
        'maximum_absolute_crosscheck_error_um': accepted_error,
        'unconverted_key_maximum_absolute_error': direct_error,
        'aligned_HE_image_spatial_used_for_distances': False,
        'fitted_affine_or_row_reordering_used': False,
        'scope': 'P08_breast_Xenium_author_preprocessing_not_a_universal_Xenium_scale',
    }
