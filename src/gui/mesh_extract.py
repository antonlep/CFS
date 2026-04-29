"""Extract solver-ready arrays from a Gmsh model."""

from dataclasses import dataclass


@dataclass
class MeshData:
    """All mesh data needed by the solver, already 0-based."""

    nodes: list[list[float]]  # [[x, y], ...]
    elements: list[list[int]]  # [[n0, n1, ..., n5], ...] 0-based
    right_edges: list[list[int]]  # [[n0, n1, n2], ...] 0-based
    left_edges: list[list[int]]  # [[n0, n1, n2], ...] 0-based
    x_fixed_nodes: list[int]  # 1-based Gmsh tags
    y_fixed_nodes: list[int]  # 1-based Gmsh tags


def get_nodes_by_name(model, group_name):
    """Return (node_tags, coordinates) for a named physical group."""
    for dim, tag in model.getPhysicalGroups():
        if model.getPhysicalName(dim, tag) == group_name:
            return model.mesh.getNodesForPhysicalGroup(dim, tag)
    raise ValueError(f"Physical group '{group_name}' not found")


def get_edge_nodes(model, group_name) -> list[list[int]]:
    """Return list of [n1, n2, n3] node tag triples for quadratic edges."""
    for dim, tag in model.getPhysicalGroups():
        if model.getPhysicalName(dim, tag) == group_name:
            entities = model.getEntitiesForPhysicalGroup(dim, tag)
            _, _, node_tags = model.mesh.getElements(1, entities[0])
            flat = node_tags[0]
            return [
                [int(flat[i]), int(flat[i + 1]), int(flat[i + 2])]
                for i in range(0, len(flat), 3)
            ]
    raise ValueError(f"Physical group '{group_name}' not found")


def get_elements_by_name(model, group_name):
    """Return flat list of element node tags for a named physical group."""
    for dim, tag in model.getPhysicalGroups():
        if model.getPhysicalName(dim, tag) == group_name:
            entities = model.getEntitiesForPhysicalGroup(dim, tag)
            all_nodes = []
            for entity in entities:
                _, _, nodes = model.mesh.getElements(dim, entity)
                for n in nodes:
                    all_nodes.extend(int(i) for i in n)
            return all_nodes
    raise ValueError(f"Physical group '{group_name}' not found")


def extract_mesh(model) -> MeshData:
    """Extract all solver-ready data from a Gmsh model."""
    # Nodes: Gmsh gives [x0, y0, z0, x1, y1, z1, ...] — reshape to [[x, y]]
    _, raw_coords = get_nodes_by_name(model, "ALL_NODES")
    nodes = [[raw_coords[i], raw_coords[i + 1]] for i in range(0, len(raw_coords), 3)]

    # Elements: reshape flat list to 6-node sublists, convert to 0-based
    raw_elems = get_elements_by_name(model, "ALL_ELEMENTS")
    elements = [
        [raw_elems[i] - 1 for i in range(j, j + 6)] for j in range(0, len(raw_elems), 6)
    ]

    # Boundary edges: convert to 0-based
    right_edges = [
        [n - 1 for n in edge] for edge in get_edge_nodes(model, "SIDE_NODES")
    ]
    left_edges = [[n - 1 for n in edge] for edge in get_edge_nodes(model, "X_FIXED")]

    # Fixed node tags (1-based — converted in boundary_conditions.py)
    x_fixed, _ = get_nodes_by_name(model, "X_FIXED")
    y_fixed, _ = get_nodes_by_name(model, "Y_FIXED")

    return MeshData(
        nodes=nodes,
        elements=elements,
        right_edges=right_edges,
        left_edges=left_edges,
        x_fixed_nodes=list(x_fixed),
        y_fixed_nodes=list(y_fixed),
    )
