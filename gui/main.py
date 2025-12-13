import sys

import pyvista as pv
import pyvistaqt as pvqt
import solver
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


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("My App")
        self.setMinimumSize(QSize(400, 400))

        layout = QVBoxLayout()
        layout1 = QVBoxLayout()
        layout2 = QVBoxLayout()
        layout.addLayout(layout1)
        layout.addLayout(layout2)

        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

        self.plotter = pvqt.QtInteractor(self)
        self.plotter.add_mesh(pv.Sphere())
        self.plotter.show_grid()
        layout1.addWidget(self.plotter.interactor)

        self.box1 = QSpinBox()
        layout2.addWidget(self.box1)

        self.box2 = QSpinBox()
        layout2.addWidget(self.box2)

        button = QPushButton("Solve")
        button.setCheckable(True)
        button.clicked.connect(self.the_button_was_clicked)
        layout2.addWidget(button)

        self.label = QLabel()
        layout2.addWidget(self.label)

    def the_button_was_clicked(self):
        result = solver.add(self.box1.value(), self.box2.value())
        self.label.setText(str(result))
        print(f"Result: {result}")


app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec()
