"""Gmsh mesh generation for a rectangular domain."""

import gmsh


def create_mesh(width: float, height: float, mesh_size: float):
    """Create a second-order triangular mesh and return the Gmsh model.

    Physical groups created:
        ALL_NODES, ALL_ELEMENTS, SIDE_NODES (right), X_FIXED (left), Y_FIXED
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("cfs_rect")

    # ---- Geometry ----
    p1 = gmsh.model.geo.addPoint(0, 0, 0, mesh_size)
    p2 = gmsh.model.geo.addPoint(width, 0, 0, mesh_size)
    p3 = gmsh.model.geo.addPoint(width, height, 0, mesh_size)
    p4 = gmsh.model.geo.addPoint(0, height, 0, mesh_size)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    s = gmsh.model.geo.addPlaneSurface([cl])

    gmsh.model.geo.synchronize()

    # ---- Mesh ----
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)

    # ---- Physical groups ----
    def add_group(dim, entities, name):
        tag = gmsh.model.addPhysicalGroup(dim, entities)
        gmsh.model.setPhysicalName(dim, tag, name)

    add_group(2, [s], "ALL_NODES")
    add_group(2, [s], "ALL_ELEMENTS")
    add_group(1, [l2], "SIDE_NODES")
    add_group(1, [l4], "X_FIXED")
    add_group(0, [p1], "Y_FIXED")

    return gmsh.model
