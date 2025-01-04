#include <pybind11/pybind11.h>
#include "../../include/compartments/Soma.hpp"

namespace py = pybind11;

void bind_Soma(py::module& m) {
  py::class_<Soma, Compartment>(m, "Soma")
    .def(py::init<const std::string&, const double&>(),
         py::arg("name")="no_name_Soma", py::arg("length")=10,
         "Create soma with certain name and length")
    .def(py::init<const std::string&, const double&, const double&, const double&, const double&, const double&>(),
    py::arg("name")="no_name_Soma", py::arg("length")=10, py::arg("x")=0, py::arg("y")=0, py::arg("z")=0, py::arg("radius")=5,
         "Create soma with certain name and length, location (x,y,z) and radius");
}
