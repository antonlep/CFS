"""PyVista grid construction and result mapping."""

import pyvista as pv


RESULT_TYPES = [
    "stress x",
    "stress y",
    "shear xy",
    "disp x",
    "disp y",
    "temperature",
]


def build_grid(nodes, elements) -> pv.UnstructuredGrid:
    """Create a PyVista unstructured grid from nodes and elements."""
    cells = []
    for e in elements:
        cells.extend([6, *e])

    celltypes = [pv.CellType.QUADRATIC_TRIANGLE] * len(elements)
    points = [[x, y, 0.0] for x, y in nodes]

    return pv.UnstructuredGrid(cells, celltypes, points)


def apply_results(grid: pv.UnstructuredGrid, results, result_name: str):
    """Map a single result field onto the grid as point data."""
    if result_name in ("stress x", "stress y", "shear xy"):
        idx = {"stress x": 0, "stress y": 1, "shear xy": 2}[result_name]
        grid.point_data[result_name] = [s[idx] for s in results.stress]

    elif result_name in ("disp x", "disp y"):
        idx = {"disp x": 0, "disp y": 1}[result_name]
        grid.point_data[result_name] = [d[idx] for d in results.disp]

    elif result_name == "temperature":
        grid.point_data[result_name] = results.temperature


def deform_grid(grid: pv.UnstructuredGrid, nodes, results, scale: float):
    """Apply scaled displacements to grid points."""
    grid.points = [
        [x + scale * ux, y + scale * uy, 0.0]
        for (x, y), (ux, uy) in zip(nodes, results.disp)
    ]
