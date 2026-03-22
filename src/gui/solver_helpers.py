import solver


def solve_from_raw(nodes, elements, fixed, fixed_values, forces):
    inp = solver.SolverInput()
    inp.nodes = [solver.Node(float(x), float(y)) for x, y in nodes]
    inp.elements = elements
    inp.fixed_dofs = fixed
    inp.fixed_values = fixed_values
    inp.forces = forces
    return solver.solve_from_data(inp)
