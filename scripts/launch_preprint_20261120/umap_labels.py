"""Direct labels for frozen categorical UMAP panels."""
from math import cos, pi, sin

import numpy as np
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.transforms import Bbox


def add_umap_labels(ax, xy, codes, labels, *, fontsize=6, nudges=None):
    """Anchor labels at category medians and place them in display space.

    Labels follow the code dictionary order. Use None to omit an entry and
    optional nudges (in points, keyed by code) to prefer a specific position.
    Displaced labels get a leader line back to their median.
    """
    nudges = nudges or {}
    figure = ax.figure
    # Resolve scatter autoscaling before measuring label positions.
    ax.get_xlim()
    ax.get_ylim()
    renderer = RendererAgg(*figure.bbox.size, figure.dpi)
    occupied = []
    annotations = []
    offsets = [(0, 0)] + [
        (radius * cos(angle * pi / 4), radius * sin(angle * pi / 4))
        for radius in range(6, 43, 4)
        for angle in range(8)
    ]

    for code, label in enumerate(labels):
        if label is None:
            continue
        points = xy[codes == code]
        if not len(points):
            continue
        center = np.median(points, axis=0)
        dx, dy = nudges.get(code, (0, 0))
        text = ax.annotate(
            label, center, xytext=(dx, dy), textcoords="offset points",
            ha="center", va="center", fontsize=fontsize, color="#202A35",
            bbox=dict(boxstyle="round,pad=.12", fc="white", ec="none", alpha=.80),
            zorder=5,
        )

        best = None
        for ox, oy in offsets:
            position = (dx + ox, dy + oy)
            text.set_position(position)
            measured = text.get_window_extent(renderer).padded(renderer.points_to_pixels(.9))
            if not all(ax.bbox.contains(*corner) for corner in (measured.min, measured.max)):
                continue
            overlaps = [Bbox.intersection(measured, other) for other in occupied]
            area = sum(box.width * box.height for box in overlaps if box is not None)
            if best is None or area < best[0]:
                best = (area, position, measured)
            if area == 0:
                break

        if best is None:
            text.set_position((0, 0))
            occupied.append(text.get_window_extent(renderer))
        else:
            _, position, measured = best
            text.set_position(position)
            occupied.append(measured)
            if position != (0, 0):
                ax.annotate(
                    "", center, xytext=position, textcoords="offset points",
                    arrowprops=dict(arrowstyle="-", color="#57616B", lw=.45),
                    zorder=4,
                )
        annotations.append(text)

    return annotations
