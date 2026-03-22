import gmsh


def write_nset(filename, name, dim, entity):
    tags, _, _ = gmsh.model.mesh.getNodes(dim, entity, includeBoundary=True)
    with open(filename, "a") as f:
        f.write(f"*NSET,NSET={name}\n")
        for i, tag in enumerate(tags):
            f.write(f"{tag}, ")
            if (i + 1) % 16 == 0:
                f.write("\n")
        f.write("\n")


def create_mesh(width, height, lc):

    filename = "mesh.inp"

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.logger.start()
    gmsh.model.add("second_order_example")

    # ---- Geometry ----
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(width, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(width, height, 0, lc)
    p4 = gmsh.model.geo.addPoint(0, height, 0, lc)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    s = gmsh.model.geo.addPlaneSurface([cl])

    gmsh.model.geo.synchronize()

    # ---- Generate 2D mesh ----
    gmsh.model.mesh.generate(2)

    # ---- Convert to second order ----
    gmsh.model.mesh.setOrder(2)

    # ---- Recommended for curved geometries ----
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)

    surf_tag = gmsh.model.addPhysicalGroup(2, [s])
    gmsh.model.setPhysicalName(2, surf_tag, "ALL_NODES")

    surf_tag = gmsh.model.addPhysicalGroup(1, [l2])
    gmsh.model.setPhysicalName(1, surf_tag, "SIDE_NODES")

    surf_tag = gmsh.model.addPhysicalGroup(0, [p1])
    gmsh.model.setPhysicalName(0, surf_tag, "Y_FIXED")

    surf_tag = gmsh.model.addPhysicalGroup(1, [l4])
    gmsh.model.setPhysicalName(1, surf_tag, "X_FIXED")

    surfaces = gmsh.model.getEntities(dim=2)
    surface_tags = [s[1] for s in surfaces]
    group = gmsh.model.addPhysicalGroup(2, surface_tags)
    gmsh.model.setPhysicalName(2, group, "ALL_ELEMENTS")

    gmsh.write(filename)

    side_nodeTags, _, _ = gmsh.model.mesh.getNodes(1, l2, includeBoundary=True)
    name = "SIDE_NODES"
    dim = 1
    entity = l2
    write_nset(filename, name, dim, entity)

    name = "Y_FIXED"
    dim = 0
    entity = p1
    write_nset(filename, name, dim, entity)

    name = "X_FIXED"
    dim = 1
    entity = l4
    write_nset(filename, name, dim, entity)

    log = gmsh.logger.get()
    gmsh.logger.stop()
    # gmsh.finalize()

    f = open("gmsh.log", "w")
    for i in log:
        f.write(f"{i}\n")
    f.close()

    return gmsh.model
