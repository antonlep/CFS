import sys

import pyvista as pv
import pyvistaqt as pvqt
from PySide6.QtCore import QSize
from PySide6.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

import solver


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("My App")
        self.setMinimumSize(QSize(600, 600))

        layout = QVBoxLayout()
        layout1 = QVBoxLayout()
        layout2 = QVBoxLayout()
        layout.addLayout(layout1)
        layout.addLayout(layout2)

        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

        label1 = QLabel("height:")
        layout2.addWidget(label1)
        self.box1 = QSpinBox()
        self.box1.setValue(10)
        layout2.addWidget(self.box1)

        label2 = QLabel("width:")
        layout2.addWidget(label2)
        self.box2 = QSpinBox()
        self.box2.setValue(20)
        layout2.addWidget(self.box2)

        label4 = QLabel("force:")
        layout2.addWidget(label4)
        self.box4 = QSpinBox()
        self.box4.setMinimum(0)
        self.box4.setMaximum(10000)
        self.box4.setValue(1000)
        self.box4.setSingleStep(1000)
        layout2.addWidget(self.box4)

        button = QPushButton("Solve")
        button.setCheckable(True)
        button.clicked.connect(self.the_button_was_clicked)
        layout2.addWidget(button)

        label5 = QLabel("displacement:")
        layout2.addWidget(label5)
        self.label = QLabel()
        layout2.addWidget(self.label)

        self.plotter = pvqt.QtInteractor(self)
        layout1.addWidget(self.plotter.interactor)

    def the_button_was_clicked(self):
        height = self.box1.value()
        width = self.box2.value()
        force = self.box4.value()
        nodes = [
            [0, 0],
            [0, height],
            [width / 2, height],
            [width / 2, 0],
            [width, 0],
            [width, height],
        ]
        elements = [[0, 2, 1], [0, 3, 2], [2, 3, 5], [3, 4, 5]]
        fixed = [0, 1, 2, 3]
        displacements = [0, 0, 0, 0]
        forces = [0, 0, 0, 0, 0, 0, 0, 0, force / 2, 0, force / 2, 0]
        cells = []
        for e in elements:
            cells.extend([3, *e])
        celltypes = [pv.CellType.TRIANGLE] * len(elements)
        points = [[x, y, 0.0] for x, y in nodes]
        grid = pv.UnstructuredGrid(cells, celltypes, points)
        results = solver.solve_from_data(nodes, elements, fixed, displacements, forces)
        cell_scalars = [i[0] for i in results]
        grid.cell_data["Stress"] = cell_scalars
        self.plotter.clear()
        self.plotter.add_axes_at_origin()
        self.plotter.add_mesh(grid, show_edges=True, scalars="Stress", cmap="rainbow")


app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec()
