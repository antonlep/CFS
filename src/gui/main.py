import sys
from dataclasses import dataclass

import pyvista as pv
import pyvistaqt as pvqt
from PySide6.QtCore import QSize, QTimer
from PySide6.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
    QComboBox,
    QFormLayout,
)

from solver_helpers import solve_from_raw


@dataclass
class GeometryParams:
    height: int
    width: int


@dataclass
class LoadParams:
    force_x: int
    force_y: int


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("CFS GUI")
        self.setMinimumSize(QSize(600, 600))

        self.build_ui()

        self.solve_timer = QTimer(self)
        self.solve_timer.setSingleShot(True)
        self.solve_timer.setInterval(200)
        self.solve_timer.timeout.connect(self.solve)

    def build_ui(self):
        root_layout = QVBoxLayout()
        controls_layout = QFormLayout()

        self.height_box = self._spinbox(default=10, step=2)
        self.width_box = self._spinbox(default=20, step=2)
        controls_layout.addRow("Height", self.height_box)
        controls_layout.addRow("Width", self.width_box)

        self.fx_box = self._spinbox(default=5000, step=5000, max_val=50000)
        self.fy_box = self._spinbox(default=5000, step=5000, max_val=50000)
        controls_layout.addRow("Force X", self.fx_box)
        controls_layout.addRow("Force Y", self.fy_box)

        self.result_combo = QComboBox()
        self.result_combo.addItems(
            ["stress x", "stress y", "shear xy", "disp x", "disp y"]
        )
        controls_layout.addRow("Plot", self.result_combo)

        self.solve_button = QPushButton("Solve")
        self.solve_button.clicked.connect(self.solve)
        controls_layout.addRow(self.solve_button)

        self.result_label = QLabel("-")
        controls_layout.addRow("Displacement", self.result_label)

        self.plotter = pvqt.QtInteractor(self)

        root_layout.addLayout(controls_layout)
        root_layout.addWidget(self.plotter.interactor)

        central = QWidget()
        central.setLayout(root_layout)
        self.setCentralWidget(central)

        self.height_box.valueChanged.connect(self.schedule_solve)
        self.width_box.valueChanged.connect(self.schedule_solve)
        self.fx_box.valueChanged.connect(self.schedule_solve)
        self.fy_box.valueChanged.connect(self.schedule_solve)
        self.result_combo.currentIndexChanged.connect(self.schedule_solve)

    @staticmethod
    def _spinbox(default=0, step=1, min_val=0, max_val=100000):
        box = QSpinBox()
        box.setRange(min_val, max_val)
        box.setSingleStep(step)
        box.setValue(default)
        return box

    def _read_geometry(self) -> GeometryParams:
        return GeometryParams(
            height=self.height_box.value(),
            width=self.width_box.value(),
        )

    def _read_loads(self) -> LoadParams:
        return LoadParams(
            force_x=self.fx_box.value(),
            force_y=self.fy_box.value(),
        )

    def solve(self):
        if self.height_box.value() <= 0 or self.width_box.value() <= 0:
            return
        geom = self._read_geometry()
        loads = self._read_loads()

        nodes, elements = self._build_mesh(geom)
        fixed, displacements, forces = self._build_bc(loads, len(nodes))

        try:
            results = solve_from_raw(
                nodes,
                elements,
                fixed,
                displacements,
                forces,
            )
        except Exception as e:
            self.result_label.setText(f"Solver error: {e}")
            return

        grid = self._build_grid(nodes, elements)
        self._apply_results(grid, results)

        self.plotter.clear()
        self.plotter.add_axes_at_origin()
        self.plotter.add_mesh(
            grid,
            show_edges=True,
            scalars=self.result_combo.currentText(),
            cmap="rainbow",
        )

    def _build_mesh(self, geom: GeometryParams):
        h, w = geom.height, geom.width

        nodes = [
            [0, 0],
            [0, h],
            [w / 2, h],
            [w / 2, 0],
            [w, 0],
            [w, h],
        ]

        elements = [
            [0, 2, 1],
            [0, 3, 2],
            [2, 3, 5],
            [3, 4, 5],
        ]

        return nodes, elements

    def _build_bc(self, loads: LoadParams, n_nodes: int):
        fixed = [0, 1, 2, 3]
        displacements = [0, 0, 0, 0]

        fx = loads.force_x / 2
        fy = loads.force_y / 2

        forces = [0.0] * (2 * n_nodes)
        forces[-4:] = [fx, fy, fx, fy]

        return fixed, displacements, forces

    def _build_grid(self, nodes, elements):
        cells = []
        for e in elements:
            cells.extend([3, *e])

        celltypes = [pv.CellType.TRIANGLE] * len(elements)
        points = [[x, y, 0.0] for x, y in nodes]

        return pv.UnstructuredGrid(cells, celltypes, points)

    def _apply_results(self, grid, results):
        name = self.result_combo.currentText()

        if name.startswith("stress"):
            idx = {"stress x": 0, "stress y": 1, "shear xy": 2}[name]
            grid.cell_data[name] = [s[idx] for s in results.stress]

        else:  # displacement
            idx = {"disp x": 0, "disp y": 1}[name]
            grid.point_data[name] = [d[idx] for d in results.disp]

    def schedule_solve(self):
        self.solve_timer.start()


# ----------------------------
# App entry point
# ----------------------------


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
