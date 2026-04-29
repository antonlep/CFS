"""Build solver boundary conditions from mesh data and load parameters."""

from dataclasses import dataclass
from solver import BoundaryEdge


@dataclass
class LoadParams:
    force_x: float
    force_y: float
    temp_left: float
    temp_right: float


@dataclass
class BoundaryConditions:
    fixed_dofs: list[int]
    fixed_values: list[float]
    forces: list[float]
    convection_bcs: list[BoundaryEdge]


def build_boundary_conditions(
    loads: LoadParams,
    n_nodes: int,
    x_fixed_tags: list[int],
    y_fixed_tags: list[int],
    right_edges: list[list[int]],
    left_edges: list[list[int]],
) -> BoundaryConditions:
    """Build all BCs from mesh data and load parameters.

    Args:
        x_fixed_tags / y_fixed_tags: 1-based Gmsh node tags
        right_edges / left_edges: 0-based node index triples
    """
    # --- Fixed DOFs (convert 1-based tags to 0-based DOF indices) ---
    fixed_dofs = [int((tag - 1) * 2) for tag in x_fixed_tags] + [
        int((tag - 1) * 2 + 1) for tag in y_fixed_tags
    ]
    fixed_values = [0.0] * len(fixed_dofs)

    # --- Distributed forces on right edge ---
    forces = [0.0] * (2 * n_nodes)
    n_edges = len(right_edges)
    if n_edges > 0:
        fx_per_edge = loads.force_x / n_edges / 6
        fy_per_edge = loads.force_y / n_edges / 6

        for edge in right_edges:
            # Simpson's rule weights: end=1, mid=4 (total=6)
            weights = {0: 1.0, 1: 4.0, 2: 1.0}
            for local_idx, global_node in enumerate(edge):
                w = weights[local_idx]
                forces[global_node * 2] += w * fx_per_edge
                forces[global_node * 2 + 1] += w * fy_per_edge

    # --- Convection BCs ---
    convection_bcs = []
    for edge in right_edges:
        convection_bcs.append(
            BoundaryEdge(
                n1=edge[0],
                n2=edge[1],
                n3=edge[2],
                h=1000.0,
                Tinf=loads.temp_right,
            )
        )
    for edge in left_edges:
        convection_bcs.append(
            BoundaryEdge(
                n1=edge[0],
                n2=edge[1],
                n3=edge[2],
                h=1000.0,
                Tinf=loads.temp_left,
            )
        )

    return BoundaryConditions(
        fixed_dofs=fixed_dofs,
        fixed_values=fixed_values,
        forces=forces,
        convection_bcs=convection_bcs,
    )
