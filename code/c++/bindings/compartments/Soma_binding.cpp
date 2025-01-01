#include <pybind11/pybind11.h>
#include "../../include/compartments/Soma.hpp"

namespace py = pybind11;

void bind_Soma(py::module& m) {
  py::class_<Soma, Compartment>(m, "Soma")
    .def(py::init<const std::string&, const double&>(),
         py::arg("name")="no_name_Soma", py::arg("length")=10,
         "Create soma with certain name and length");
}
