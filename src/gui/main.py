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

        label3 = QLabel("thickness:")
        layout2.addWidget(label3)
        self.box3 = QSpinBox()
        self.box3.setValue(1)
        layout2.addWidget(self.box3)

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
        box4 = self.box4.value()
        nodes = [[0, 0], [0, 10], [20, 10], [20, 0], [30, 0], [30, 10]]
        elements = [[0, 2, 1], [0, 3, 2], [2, 3, 5], [3, 4, 5]]
        fixed = [0, 1, 2, 3]
        displacements = [0, 0, 0, 0]
        forces = [0, 0, 0, 0, 0, 0, 0, 0, box4, 0, box4, 0]
        cells = []
        for e in elements:
            cells.extend([3, *e])
        celltypes = [pv.CellType.TRIANGLE] * len(elements)
        points = [[x, y, 0.0] for x, y in nodes]
        grid = pv.UnstructuredGrid(cells, celltypes, points)
        # grid.plot(show_edges=True)
        result = solver.solve_from_data(nodes, elements, fixed, displacements, forces)
        self.label.setText(str(result))
        print(f"Result: {result}")
        # grid = pv.ImageData()
        # grid.dimensions = (box1 + 1, box2 + 1, box3 + 1)

        # grid.origin = (0, 0, 0)  # The bottom left corner of the data set
        # grid.spacing = (1, 1, 1)  # These are the cell sizes along each axis
        self.plotter.clear()
        self.plotter.add_axes_at_origin()
        self.plotter.add_mesh(grid, show_edges=True)


app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec()
