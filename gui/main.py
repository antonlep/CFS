import sys

import solver

# from gui import show_solution

# mesh = load_mesh("triangle.msh")
# u = solver.solve_poisson(mesh)
# show_solution(mesh, u)

u = solver.add(2, 3)
print(u)
