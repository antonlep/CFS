import sys

import numpy as np
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

        self.box1 = QSpinBox()
        self.box1.setValue(4)
        layout2.addWidget(self.box1)

        self.box2 = QSpinBox()
        self.box2.setValue(4)
        layout2.addWidget(self.box2)

        button = QPushButton("Solve")
        button.setCheckable(True)
        button.clicked.connect(self.the_button_was_clicked)
        layout2.addWidget(button)

        self.label = QLabel()
        layout2.addWidget(self.label)

        self.plotter = pvqt.QtInteractor(self)
        layout1.addWidget(self.plotter.interactor)

    def the_button_was_clicked(self):
        box1 = self.box1.value()
        box2 = self.box2.value()
        result = solver.add(box1, box2)
        self.label.setText(str(result))
        print(f"Result: {result}")
        # Create the spatial reference
        grid = pv.ImageData()

        values = np.linspace(0, 1, box1 * box2).reshape((box1, box2, 1))
        values.shape

        # Set the grid dimensions: shape + 1 because we want to inject our values on
        #   the CELL data
        grid.dimensions = np.array(values.shape) + 1

        # Edit the spatial reference
        grid.origin = (0, 0, 0)  # The bottom left corner of the data set
        grid.spacing = (1, 1, 1)  # These are the cell sizes along each axis

        # Now plot the grid
        # grid.plot(show_edges=True)
        self.plotter.clear()
        # self.plotter.show_grid()
        self.plotter.add_mesh(grid, show_edges=True)
        # Add the data values to the cell data
        # grid.cell_data["values"] = values.flatten(order="F")  # Flatten the array


app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec()
