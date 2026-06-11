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
    QHBoxLayout,
    QWidget,
    QDoubleSpinBox,
    QTabWidget,
    QGroupBox,
)

from boundary_conditions import (
    LoadParams,
    MaterialParams,
    build_boundary_conditions,
    build_material,
)
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

    def _build_ui(self):
        """Build the main window UI."""
        # Build all input forms
        geom_form = self._build_geometry_form()
        mat_form = self._build_material_form()
        load_form = self._build_loads_form()
        display_form = self._build_display_form()
        output_form = self._build_output_form()

        # Assemble layout
        central = self._assemble_layout(
            geom_form, mat_form, load_form, display_form, output_form
        )
        self.setCentralWidget(central)

        # Wire up signals
        self._connect_signals()

    # ── Form Builders ──────────────────────────────────────────────

    def _build_geometry_form(self) -> QGroupBox:
        self._height_box = self._spinbox(default=10, step=2)
        self._width_box = self._spinbox(default=20, step=2)
        self._t_box = self._double_spinbox(default=1, step=0.1)

        form = QFormLayout()
        form.addRow("Height", self._height_box)
        form.addRow("Width", self._width_box)
        form.addRow("Thickness", self._t_box)

        box = QGroupBox("Geometry")
        box.setLayout(form)
        return box

    def _build_material_form(self) -> QGroupBox:
        self._e_box = self._spinbox(default=210, step=10)
        self._nu_box = self._double_spinbox(default=0.3, step=0.05)
        self._k_box = self._spinbox(default=45, step=1)
        self._alpha_box = self._double_spinbox(default=12, step=1)

        form = QFormLayout()
        form.addRow("Elastic modulus", self._e_box)
        form.addRow("Poisson's ratio", self._nu_box)
        form.addRow("Thermal conductivity", self._k_box)
        form.addRow("Thermal expansion coeff", self._alpha_box)

        box = QGroupBox("Material")
        box.setLayout(form)
        return box

    def _build_loads_form(self) -> QGroupBox:
        self._fx_box = self._spinbox(default=10, step=10, max_val=10_000)
        self._fy_box = self._spinbox(default=10, step=10, max_val=10_000)
        self._temp_left_box = self._spinbox(
            default=300, step=10, min_val=200, max_val=400
        )
        self._temp_right_box = self._spinbox(
            default=300, step=10, min_val=200, max_val=400
        )
        self._t0_box = self._spinbox(default=273, step=10)

        form = QFormLayout()
        form.addRow("Distr load X", self._fx_box)
        form.addRow("Distr load Y", self._fy_box)
        form.addRow("Temperature left", self._temp_left_box)
        form.addRow("Temperature right", self._temp_right_box)
        form.addRow("Reference temperature", self._t0_box)

        box = QGroupBox("Loads")
        box.setLayout(form)
        return box

    def _build_display_form(self) -> QGroupBox:
        self._mesh_size_box = self._double_spinbox(
            default=5, step=0.1, min_val=0.1, max_val=100, decimals=1
        )
        self._scale_box = self._spinbox(default=100, step=50, max_val=10_000)
        self._max_box = self._spinbox(
            default=100, step=10, min_val=-10_000, max_val=10_000
        )
        self._min_box = self._spinbox(
            default=-100, step=10, min_val=-10_000, max_val=10_000
        )

        self._auto_scale_box = QCheckBox()
        self._auto_scale_box.setChecked(True)
        self._show_mesh_box = QCheckBox()

        self._result_combo = QComboBox()
        self._result_combo.addItems(RESULT_TYPES)

        form = QFormLayout()
        form.addRow("Plot", self._result_combo)
        form.addRow("Auto scale", self._auto_scale_box)
        form.addRow("Scale max", self._max_box)
        form.addRow("Scale min", self._min_box)
        form.addRow("Mesh size", self._mesh_size_box)
        form.addRow("Display mesh", self._show_mesh_box)
        form.addRow("Deformation scale", self._scale_box)

        box = QGroupBox("Display")
        box.setLayout(form)
        return box

    def _build_output_form(self) -> QFormLayout:
        solve_btn = QPushButton("Solve")
        solve_btn.clicked.connect(self._solve)

        self._result_label = QLabel("-")
        self._plotter = pvqt.QtInteractor(self)

        form = QFormLayout()
        form.addRow(solve_btn)
        form.addRow("Status", self._result_label)
        form.addRow(self._plotter.interactor)

        return form

    # ── Layout Assembly ────────────────────────────────────────────

    def _assemble_layout(
        self,
        geom_form: QFormLayout,
        mat_form: QFormLayout,
        load_form: QFormLayout,
        display_form: QFormLayout,
        output_form: QFormLayout,
    ) -> QWidget:
        """Compose all forms into the main window layout."""
        # Tab 1: interactive inputs
        interactive_layout = QHBoxLayout()
        interactive_layout.addWidget(geom_form)
        interactive_layout.addWidget(mat_form)
        interactive_layout.addWidget(load_form)

        interactive_tab = QWidget()
        interactive_tab.setLayout(interactive_layout)

        # Tab 2: solve input (placeholder)
        solve_input_tab = QWidget()
        solve_input_tab.setLayout(QHBoxLayout())

        tabs = QTabWidget()
        tabs.addTab(interactive_tab, "interactive")
        tabs.addTab(solve_input_tab, "solve input")

        # Top: tabs on the left, display settings on the right
        top_row = QHBoxLayout()
        top_row.addWidget(tabs)
        top_row.addWidget(display_form)

        # Full window: top row + output below
        root = QVBoxLayout()
        root.addLayout(top_row)
        root.addLayout(output_form)

        central = QWidget()
        central.setLayout(root)
        return central

    # ── Signal Wiring ──────────────────────────────────────────────

    def _connect_signals(self):
        """Connect all input widgets to the debounced solve trigger."""
        value_widgets = [
            self._height_box,
            self._width_box,
            self._t_box,
            self._e_box,
            self._nu_box,
            self._k_box,
            self._alpha_box,
            self._fx_box,
            self._fy_box,
            self._temp_left_box,
            self._temp_right_box,
            self._t0_box,
            self._mesh_size_box,
            self._scale_box,
            self._max_box,
            self._min_box,
        ]
        for w in value_widgets:
            w.valueChanged.connect(self._schedule_solve)

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

    def _read_material(self) -> MaterialParams:
        return MaterialParams(
            e=self._e_box.value(),
            nu=self._nu_box.value(),
            t=self._t_box.value(),
            k=self._k_box.value(),
            alpha=self._alpha_box.value(),
            t0=self._t0_box.value(),
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
        material_params = self._read_material()

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

        material = build_material(material_params)

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
                material,
            )
        except Exception as e:
            self._result_label.setText(f"Solver error: {e}")
            return

        self._result_label.setText("OK")

        # 4. Visualize
        result_name = self._result_combo.currentText()
        grid = build_grid(mesh.nodes, mesh.elements, material_params.t)
        deform_grid(grid, mesh.nodes, results, display.scale, material_params.t)
        apply_results(grid, results, result_name)

        print(f"Number of elements: {len(mesh.elements)}")
        print(f"Number of nodes: {len(mesh.nodes)}")

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
                feature_edges=True,
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

    @staticmethod
    def _double_spinbox(default=0, step=1, min_val=0, max_val=100000, decimals=2):
        box = QDoubleSpinBox()
        box.setRange(min_val, max_val)
        box.setSingleStep(step)
        box.setValue(default)
        box.setDecimals(decimals)
        return box


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
