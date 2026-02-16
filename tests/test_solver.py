from src.gui.solver_helpers import solve_from_raw


def test_two_element():
    nodes = [
        [0, 0],
        [0, 5],
        [0, 10],
        [5, 0],
        [5, 5],
        [5, 10],
        [10, 0],
        [10, 5],
        [10, 10],
        [15, 0],
        [15, 5],
        [15, 10],
        [20, 0],
        [20, 5],
        [20, 10],
    ]
    elements = [
        [0, 8, 2, 4, 5, 1],
        [0, 6, 8, 3, 7, 4],
        [6, 12, 14, 9, 13, 10],
        [6, 14, 8, 10, 11, 7],
    ]
    fixed = [0, 1, 2, 4]
    displacements = [0, 0, 0, 0]
    forces = [0.0] * (2 * len(nodes))
    forces[-6:] = [1666, 0, 6666, 0, 1666, 0]
    result = solve_from_raw(nodes, elements, fixed, displacements, forces)

    print(result.stress)

    for n in result.stress:
        assert n[0] > 999
        assert n[0] < 1001
        assert abs(n[1]) < 0.2
        assert abs(n[2]) < 0.2

    assert result.disp[12][0] > 9.5e-5
    assert result.disp[13][0] > 9.5e-5
    assert result.disp[14][0] > 9.5e-5
    assert result.disp[12][0] < 9.525e-5
    assert result.disp[13][0] < 9.525e-5
    assert result.disp[14][0] < 9.525e-5
    for n in result.disp:
        assert abs(n[1]) < 2e-5
        assert abs(n[1]) < 2e-5
        assert abs(n[1]) < 2e-5
