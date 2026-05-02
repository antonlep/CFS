import solver


def solve_from_raw(
    nodes, elements, fixed, fixed_values, forces, convection_bcs, traction_bcs
):
    inp = solver.SolverInput()
    inp.nodes = [solver.Node(float(x), float(y)) for x, y in nodes]
    inp.elements = elements
    inp.fixed_dofs = fixed
    inp.fixed_values = fixed_values
    inp.forces = forces
    inp.convection_bcs = convection_bcs
    inp.traction_bcs = traction_bcs
    E = 210e6
    nu = 0.3
    t = 1.0
    k = 45.0
    alpha = 12e-6
    T0 = 273.0
    inp.material = solver.MaterialProperties(E=E, nu=nu, t=t, k=k, alpha=alpha, T0=T0)
    return solver.solve_from_data(inp)
