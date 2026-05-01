#include "solver_api.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_solver /*unused*/, m /*unused*/) {
  py::class_<Node>(m, "Node")
      .def(py::init<double, double>())
      .def_readwrite("x", &Node::x)
      .def_readwrite("y", &Node::y);

  py::class_<BoundaryEdge>(m, "BoundaryEdge")
      .def(py::init<>())
      .def(py::init([](int n1, int n2, int n3, double h, double Tinf) {
             BoundaryEdge e;
             e.n1 = n1;
             e.n2 = n2;
             e.n3 = n3;
             e.h = h;
             e.Tinf = Tinf;
             return e;
           }),
           py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("h"),
           py::arg("Tinf"))
      .def_readwrite("n1", &BoundaryEdge::n1)
      .def_readwrite("n2", &BoundaryEdge::n2)
      .def_readwrite("n3", &BoundaryEdge::n3)
      .def_readwrite("h", &BoundaryEdge::h)
      .def_readwrite("Tinf", &BoundaryEdge::Tinf);

  py::class_<TractionEdge>(m, "TractionEdge")
      .def(py::init<>())
      .def(py::init([](int n1, int n2, int n3, double tx, double ty) {
             return TractionEdge{n1, n2, n3, tx, ty};
           }),
           py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("tx"),
           py::arg("ty"))
      .def_readwrite("n1", &TractionEdge::n1)
      .def_readwrite("n2", &TractionEdge::n2)
      .def_readwrite("n3", &TractionEdge::n3)
      .def_readwrite("tx", &TractionEdge::tx)
      .def_readwrite("ty", &TractionEdge::ty);

  py::class_<SolverInput>(m, "SolverInput")
      .def(py::init<>())
      .def_readwrite("nodes", &SolverInput::nodes)
      .def_readwrite("elements", &SolverInput::elements)
      .def_readwrite("fixed_dofs", &SolverInput::fixed_dofs)
      .def_readwrite("fixed_values", &SolverInput::fixed_values)
      .def_readwrite("forces", &SolverInput::forces)
      .def_readwrite("convection_bcs", &SolverInput::boundary_edges)
      .def_readwrite("traction_bcs", &SolverInput::traction_edges);

  py::class_<SolverOutput>(m, "SolverOutput")
      .def_readonly("stress", &SolverOutput::stress)
      .def_readonly("disp", &SolverOutput::disp)
      .def_readonly("temperature", &SolverOutput::temperature);

  m.def("solve_from_data", &solve_from_data);

  m.def("solve_from_file", &solve_from_file);
}
