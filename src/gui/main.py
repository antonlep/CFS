import sys
from dataclasses import dataclass
from mesh import create_mesh
from solver import BoundaryEdge

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
    QCheckBox,
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
    temp_left: int
    temp_right: int


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

        self.fx_box = self._spinbox(default=5000, step=5000, max_val=1000000)
        self.fy_box = self._spinbox(default=5000, step=5000, max_val=1000000)
        self.temp_left_box = self._spinbox(
            default=300, step=1, min_val=280, max_val=320
        )
        self.temp_right_box = self._spinbox(
            default=300, step=1, min_val=280, max_val=320
        )
        self.size_box = self._spinbox(default=10, step=1, max_val=100)
        self.max_box = self._spinbox(
            default=100, step=10, max_val=10000, min_val=-10000
        )
        self.min_box = self._spinbox(
            default=-100, step=10, max_val=10000, min_val=-10000
        )
        controls_layout.addRow("Force X", self.fx_box)
        controls_layout.addRow("Force Y", self.fy_box)
        controls_layout.addRow("Temperature left", self.temp_left_box)
        controls_layout.addRow("Temperature right", self.temp_right_box)
        controls_layout.addRow("Mesh size", self.size_box)

        self.result_combo = QComboBox()
        self.result_combo.addItems(
            ["stress x", "stress y", "shear xy", "disp x", "disp y", "temperature"]
        )
        self.scale_box = self._spinbox(default=1000, step=500, max_val=10000)
        self.auto_box = QCheckBox()
        self.auto_box.setChecked(True)
        self.mesh_box = QCheckBox()
        controls_layout.addRow("Scale factor", self.scale_box)
        controls_layout.addRow("Scale max", self.max_box)
        controls_layout.addRow("Scale min", self.min_box)
        controls_layout.addRow("Auto scale", self.auto_box)
        controls_layout.addRow("Display mesh", self.mesh_box)
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
        self.temp_left_box.valueChanged.connect(self.schedule_solve)
        self.temp_right_box.valueChanged.connect(self.schedule_solve)
        self.size_box.valueChanged.connect(self.schedule_solve)
        self.result_combo.currentIndexChanged.connect(self.schedule_solve)
        self.scale_box.valueChanged.connect(self.schedule_solve)
        self.mesh_box.stateChanged.connect(self.schedule_solve)
        self.auto_box.stateChanged.connect(self.schedule_solve)
        self.max_box.valueChanged.connect(self.schedule_solve)
        self.min_box.valueChanged.connect(self.schedule_solve)

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
            temp_left=self.temp_left_box.value(),
            temp_right=self.temp_right_box.value(),
        )

    def _read_options(self) -> LoadParams:
        return (
            self.scale_box.value(),
            self.mesh_box.isChecked(),
            self.size_box.value(),
            self.min_box.value(),
            self.max_box.value(),
            self.auto_box.isChecked(),
        )

    def _get_nodes_by_name(self, model, group_name):
        # Get all physical groups: list of (dim, tag)
        groups = model.getPhysicalGroups()

        for dim, tag in groups:
            # Check if this tag's name matches our target
            if model.getPhysicalName(dim, tag) == group_name:
                # Returns (nodeTags, coord, parametricCoord)
                nodes = model.mesh.getNodesForPhysicalGroup(dim, tag)
                return nodes

    def _get_edge_nodes(self, model, group_name):
        groups = model.getPhysicalGroups()
        nodes = []
        coords = []
        for dim, tag in groups:
            if model.getPhysicalName(dim, tag) == group_name:
                entities = model.getEntitiesForPhysicalGroup(dim, tag)
                types, elementTags, nodeTags = model.mesh.getElements(1, entities[0])
                nodeTags = nodeTags[0]
                nodeTags = [
                    nodeTags[i : i + 3] for i in range(0, len(nodeTags), 3)  # noqa
                ]
                for n in nodeTags:
                    node_coords = []
                    node_nodes = []
                    for i in n:
                        c = model.mesh.getNode(i)
                        node_coords.append(float(c[0][1]))
                        node_nodes.append(int(i))
                    coords.append(node_coords)
                    nodes.append(node_nodes)
        return elementTags, nodes, coords

    def _get_elements_by_name(self, model, group_name):
        groups = model.getPhysicalGroups()

        for dim, tag in groups:
            if model.getPhysicalName(dim, tag) == group_name:

                entities = model.getEntitiesForPhysicalGroup(dim, tag)

                element_tags = []
                element_nodes = []

                for entity in entities:
                    types, tags, nodes = model.mesh.getElements(dim, entity)

                    for t in tags:
                        tt = [int(i) for i in t]
                        element_tags.extend(tt)

                    for n in nodes:
                        nn = [int(i) for i in n]
                        element_nodes.extend(nn)

                return element_tags, element_nodes

    def solve(self):
        if self.height_box.value() <= 0 or self.width_box.value() <= 0:
            return
        geom = self._read_geometry()
        loads = self._read_loads()
        scale, show_mesh, mesh_size, scale_min, scale_max, auto_scale = (
            self._read_options()
        )

        nodes, elements = self._build_mesh(geom)
        model = create_mesh(geom.width, geom.height, mesh_size)
        nn, nodes = self._get_nodes_by_name(model, "ALL_NODES")
        ee, elements = self._get_elements_by_name(model, "ALL_ELEMENTS")
        element_tags, node_tags, node_coords = self._get_edge_nodes(model, "SIDE_NODES")
        left_element_tags, left_node_tags, left_node_coords = self._get_edge_nodes(
            model, "X_FIXED"
        )
        print(element_tags)
        print(node_tags)
        print(node_coords)
        print(left_element_tags)
        print(left_node_tags)
        print(left_node_coords)
        y_fixed, _ = self._get_nodes_by_name(model, "Y_FIXED")
        x_fixed, _ = self._get_nodes_by_name(model, "X_FIXED")

        nodes = [[nodes[i], nodes[i + 1]] for i in range(0, len(nodes), 3)]
        elements = [elements[i : i + 6] for i in range(0, len(elements), 6)]  # noqa
        elements = [[n - 1 for n in e] for e in elements]

        fixed, displacements, forces, convection_bcs = self._build_bc(
            loads,
            node_tags,
            node_coords,
            left_node_tags,
            x_fixed,
            y_fixed,
            len(nodes),
        )

        print("nodes:")
        print(len(nodes))
        print("elements:")
        print(len(elements))
        print()
        try:
            results = solve_from_raw(
                nodes, elements, fixed, displacements, forces, convection_bcs
            )
        except Exception as e:
            self.result_label.setText(f"Solver error: {e}")
            return

        grid = self._build_grid(nodes, elements)
        points = []
        for (x, y), (ux, uy) in zip(nodes, results.disp):
            points.append([x + scale * ux, y + scale * uy, 0.0])

        grid.points = points
        self._apply_results(grid, results)

        self.plotter.clear()
        self.plotter.add_axes_at_origin()
        if auto_scale:
            self.plotter.add_mesh(
                grid,
                show_edges=False,
                scalars=self.result_combo.currentText(),
                cmap="rainbow",
            )
        else:
            self.plotter.add_mesh(
                grid,
                show_edges=False,
                scalars=self.result_combo.currentText(),
                cmap="rainbow",
                clim=[scale_min, scale_max],
            )

        element_edges = grid.extract_all_edges()
        outer_edges = grid.extract_feature_edges(
            boundary_edges=True,
            feature_edges=False,
            manifold_edges=False,
            non_manifold_edges=False,
        )
        if show_mesh:
            self.plotter.add_mesh(element_edges, color="black", line_width=2)
        else:
            self.plotter.add_mesh(outer_edges, color="black", line_width=2)

    def _build_mesh(self, geom: GeometryParams):
        h, w = geom.height, geom.width

        nodes = [
            [0, 0],
            [0, h / 2],
            [0, h],
            [w / 4, 0],
            [w / 4, h / 2],
            [w / 4, h],
            [2 * w / 4, 0],
            [2 * w / 4, h / 2],
            [2 * w / 4, h],
            [3 * w / 4, 0],
            [3 * w / 4, h / 2],
            [3 * w / 4, h],
            [4 * w / 4, 0],
            [4 * w / 4, h / 2],
            [4 * w / 4, h],
        ]
        elements = [
            [0, 8, 2, 4, 5, 1],
            [0, 6, 8, 3, 7, 4],
            [6, 12, 14, 9, 13, 10],
            [6, 14, 8, 10, 11, 7],
        ]

        return nodes, elements

    def _build_bc(
        self,
        loads: LoadParams,
        node_tags,
        node_coords,
        left_node_tags,
        x_fixed,
        y_fixed,
        n_nodes: int,
    ):
        x_fixed_trans = [int((i - 1) * 2) for i in x_fixed]
        y_fixed_trans = [int((i - 1) * 2 + 1) for i in y_fixed]
        fixed = x_fixed_trans + y_fixed_trans
        displacements = [0] * len(fixed)
        min_y = 1e99
        max_y = -1e99
        for e in node_coords:
            for n in e:
                if n < min_y:
                    min_y = n
                if n > max_y:
                    max_y = n
        # height = max_y - min_y
        forces = [0.0] * (2 * n_nodes)
        for i in range(0, len(node_coords)):
            local_coords = node_coords[i]
            min_node = local_coords.index(min(local_coords))
            max_node = local_coords.index(max(local_coords))
            # dist = local_coords[max_node] - local_coords[min_node]
            # fx = loads.force_x * dist / height / 6
            # fy = loads.force_y * dist / height / 6
            fx = loads.force_x / len(node_coords) / 6
            fy = loads.force_y / len(node_coords) / 6
            for c in [0, 1, 2]:
                label = node_tags[i][c]
                if c == min_node or c == max_node:
                    forces[(label - 1) * 2] += fx
                    forces[(label - 1) * 2 + 1] += fy
                else:
                    forces[(label - 1) * 2] += 4 * fx
                    forces[(label - 1) * 2 + 1] += 4 * fy
        convection_bcs = []
        for i in node_tags:
            convection_bcs.append(
                BoundaryEdge(
                    n1=int(i[0] - 1),
                    n2=int(i[1] - 1),
                    n3=int(i[2] - 1),
                    h=1000.0,
                    Tinf=float(loads.temp_right),
                )
            )
        for i in left_node_tags:
            convection_bcs.append(
                BoundaryEdge(
                    n1=int(i[0] - 1),
                    n2=int(i[1] - 1),
                    n3=int(i[2] - 1),
                    h=1000.0,
                    Tinf=float(loads.temp_left),
                )
            )
        return fixed, displacements, forces, convection_bcs

    def _build_grid(self, nodes, elements):
        cells = []
        for e in elements:
            cells.extend([6, *e])

        celltypes = [pv.CellType.QUADRATIC_TRIANGLE] * len(elements)
        points = [[x, y, 0.0] for x, y in nodes]

        return pv.UnstructuredGrid(cells, celltypes, points)

    def _apply_results(self, grid, results):
        name = self.result_combo.currentText()

        if name.startswith("stress") or name.startswith("shear"):
            idx = {"stress x": 0, "stress y": 1, "shear xy": 2}[name]
            grid.point_data[name] = [s[idx] for s in results.stress]

        elif name.startswith("disp"):  # displacement
            idx = {"disp x": 0, "disp y": 1}[name]
            grid.point_data[name] = [d[idx] for d in results.disp]
        else:
            grid.point_data[name] = results.temperature

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
