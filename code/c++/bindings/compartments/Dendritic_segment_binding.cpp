#include <pybind11/pybind11.h>
#include "../../include/compartments/Dendritic_segment.hpp"

namespace py = pybind11;

void bind_Dendritic_segment(py::module& m) {
  py::class_<Dendritic_segment, Compartment>(m, "Dendritic_segment")
    .def(py::init<Compartment&, const std::string&, const double&>());
}
