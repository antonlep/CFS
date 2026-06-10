"""PyVista grid construction and result mapping."""

import pyvista as pv
import numpy as np

RESULT_TYPES = [
    "stress x",
    "stress y",
    "shear xy",
    "disp x",
    "disp y",
    "temperature",
]


def build_grid(nodes, elements, thickness: float) -> pv.UnstructuredGrid:
    """Create 3D grid by extruding quadratic triangles into quadratic wedges."""
    n = len(nodes)
    h = thickness / 2.0

    # Three layers: bottom (z=-h), top (z=+h), middle (z=0)
    points = (
        [[x, y, -h] for x, y in nodes]  # 0      .. n-1    bottom
        + [[x, y, h] for x, y in nodes]  # n      .. 2n-1   top
        + [[x, y, 0.0] for x, y in nodes]  # 2n     .. 3n-1   mid (vertical edges)
    )

    cells = []
    celltypes = [pv.CellType.QUADRATIC_WEDGE] * len(elements)
    for e in elements:
        cells.append(15)
        cells.extend(
            [
                e[0],
                e[1],
                e[2],  # bottom corners
                e[0] + n,
                e[1] + n,
                e[2] + n,  # top corners
                e[3],
                e[4],
                e[5],  # bottom mid-edges
                e[3] + n,
                e[4] + n,
                e[5] + n,  # top mid-edges
                e[0] + 2 * n,
                e[1] + 2 * n,
                e[2] + 2 * n,  # vertical mid-edges
            ]
        )

    return pv.UnstructuredGrid(cells, celltypes, points)


def apply_results(grid: pv.UnstructuredGrid, results, result_name: str):
    """Map result values onto grid points (replicated across 3 z-layers)."""
    if result_name in ("stress x", "stress y", "shear xy"):
        idx = {"stress x": 0, "stress y": 1, "shear xy": 2}[result_name]
        values = [s[idx] for s in results.stress]
    elif result_name in ("disp x", "disp y"):
        idx = {"disp x": 0, "disp y": 1}[result_name]
        values = [d[idx] for d in results.disp]
    elif result_name == "temperature":
        values = list(results.temperature)

    grid.point_data[result_name] = values * 3


def deform_grid(
    grid: pv.UnstructuredGrid, nodes, results, scale: float, thickness: float
):
    """Apply scaled displacements while preserving z-layers."""
    h = thickness / 2.0
    z_levels = [-h, h, 0.0]  # bottom, top, middle

    new_points = []
    for z in z_levels:
        for (x, y), (ux, uy) in zip(nodes, results.disp):
            new_points.append([x + scale * ux, y + scale * uy, z])

    grid.points = np.array(new_points)
