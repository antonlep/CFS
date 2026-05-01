"""CFS GUI — main window and application entry point."""

import sys
from dataclasses import dataclass

import pyvistaqt as pvqt
from PySide6.QtCore import QSize, QTimer
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFormLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from boundary_conditions import LoadParams, build_boundary_conditions
from mesh import create_mesh
from mesh_extract import extract_mesh
from solver_helpers import solve_from_raw
from visualization import RESULT_TYPES, apply_results, build_grid, deform_grid


@dataclass
class GeometryParams:
    height: int
    width: int


@dataclass
class DisplayOptions:
    scale: int
    show_mesh: bool
    mesh_size: int
    scale_min: int
    scale_max: int
    auto_scale: bool


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CFS GUI")
        self.setMinimumSize(QSize(600, 600))
        self._build_ui()

        self._solve_timer = QTimer(self)
        self._solve_timer.setSingleShot(True)
        self._solve_timer.setInterval(200)
        self._solve_timer.timeout.connect(self._solve)

    # ── UI Construction ──

    def _build_ui(self):
        form = QFormLayout()

        self._height_box = self._spinbox(default=10, step=2)
        self._width_box = self._spinbox(default=20, step=2)
        form.addRow("Height", self._height_box)
        form.addRow("Width", self._width_box)

        self._fx_box = self._spinbox(default=5000, step=5000, max_val=1000000)
        self._fy_box = self._spinbox(default=5000, step=5000, max_val=1000000)
        self._temp_left_box = self._spinbox(
            default=300, step=1, min_val=280, max_val=320
        )
        self._temp_right_box = self._spinbox(
            default=300, step=1, min_val=280, max_val=320
        )
        self._mesh_size_box = self._spinbox(default=10, step=1, max_val=100)
        form.addRow("Force X", self._fx_box)
        form.addRow("Force Y", self._fy_box)
        form.addRow("Temperature left", self._temp_left_box)
        form.addRow("Temperature right", self._temp_right_box)
        form.addRow("Mesh size", self._mesh_size_box)

        self._scale_box = self._spinbox(default=1000, step=500, max_val=10000)
        self._max_box = self._spinbox(
            default=100, step=10, max_val=10000, min_val=-10000
        )
        self._min_box = self._spinbox(
            default=-100, step=10, max_val=10000, min_val=-10000
        )
        self._auto_scale_box = QCheckBox()
        self._auto_scale_box.setChecked(True)
        self._show_mesh_box = QCheckBox()
        form.addRow("Scale factor", self._scale_box)
        form.addRow("Scale max", self._max_box)
        form.addRow("Scale min", self._min_box)
        form.addRow("Auto scale", self._auto_scale_box)
        form.addRow("Display mesh", self._show_mesh_box)

        self._result_combo = QComboBox()
        self._result_combo.addItems(RESULT_TYPES)
        form.addRow("Plot", self._result_combo)

        solve_btn = QPushButton("Solve")
        solve_btn.clicked.connect(self._solve)
        form.addRow(solve_btn)

        self._result_label = QLabel("-")
        form.addRow("Status", self._result_label)

        self._plotter = pvqt.QtInteractor(self)

        root = QVBoxLayout()
        root.addLayout(form)
        root.addWidget(self._plotter.interactor)

        central = QWidget()
        central.setLayout(root)
        self.setCentralWidget(central)

        # Connect all controls to debounced solve
        for widget in [
            self._height_box,
            self._width_box,
            self._fx_box,
            self._fy_box,
            self._temp_left_box,
            self._temp_right_box,
            self._mesh_size_box,
            self._scale_box,
            self._max_box,
            self._min_box,
        ]:
            widget.valueChanged.connect(self._schedule_solve)

        self._result_combo.currentIndexChanged.connect(self._schedule_solve)
        self._show_mesh_box.stateChanged.connect(self._schedule_solve)
        self._auto_scale_box.stateChanged.connect(self._schedule_solve)

    # ── Read UI State ──

    def _read_geometry(self) -> GeometryParams:
        return GeometryParams(
            height=self._height_box.value(),
            width=self._width_box.value(),
        )

    def _read_loads(self) -> LoadParams:
        return LoadParams(
            force_x=self._fx_box.value(),
            force_y=self._fy_box.value(),
            temp_left=self._temp_left_box.value(),
            temp_right=self._temp_right_box.value(),
        )

    def _read_display(self) -> DisplayOptions:
        return DisplayOptions(
            scale=self._scale_box.value(),
            show_mesh=self._show_mesh_box.isChecked(),
            mesh_size=self._mesh_size_box.value(),
            scale_min=self._min_box.value(),
            scale_max=self._max_box.value(),
            auto_scale=self._auto_scale_box.isChecked(),
        )

    # ── Solve + Plot ──

    def _solve(self):
        geom = self._read_geometry()
        if geom.height <= 0 or geom.width <= 0:
            return

        loads = self._read_loads()
        display = self._read_display()

        # 1. Mesh
        model = create_mesh(geom.width, geom.height, display.mesh_size)
        mesh = extract_mesh(model)

        # 2. Boundary conditions
        bcs = build_boundary_conditions(
            loads=loads,
            n_nodes=len(mesh.nodes),
            x_fixed_tags=mesh.x_fixed_nodes,
            y_fixed_tags=mesh.y_fixed_nodes,
            right_edges=mesh.right_edges,
            left_edges=mesh.left_edges,
        )

        # 3. Solve
        try:
            results = solve_from_raw(
                mesh.nodes,
                mesh.elements,
                bcs.fixed_dofs,
                bcs.fixed_values,
                bcs.forces,
                bcs.convection_bcs,
                bcs.traction_bcs,
            )
        except Exception as e:
            self._result_label.setText(f"Solver error: {e}")
            return

        self._result_label.setText("OK")

        # 4. Visualize
        result_name = self._result_combo.currentText()
        grid = build_grid(mesh.nodes, mesh.elements)
        deform_grid(grid, mesh.nodes, results, display.scale)
        apply_results(grid, results, result_name)

        self._plot(grid, result_name, display)

    def _plot(self, grid, result_name: str, display: DisplayOptions):
        self._plotter.clear()
        self._plotter.add_axes_at_origin()

        clim = None if display.auto_scale else [display.scale_min, display.scale_max]
        self._plotter.add_mesh(
            grid,
            show_edges=False,
            scalars=result_name,
            cmap="rainbow",
            clim=clim,
        )

        if display.show_mesh:
            edges = grid.extract_all_edges()
        else:
            edges = grid.extract_feature_edges(
                boundary_edges=True,
                feature_edges=False,
                manifold_edges=False,
                non_manifold_edges=False,
            )
        self._plotter.add_mesh(edges, color="black", line_width=2)

    def _schedule_solve(self):
        self._solve_timer.start()

    @staticmethod
    def _spinbox(default=0, step=1, min_val=0, max_val=100000):
        box = QSpinBox()
        box.setRange(min_val, max_val)
        box.setSingleStep(step)
        box.setValue(default)
        return box


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
