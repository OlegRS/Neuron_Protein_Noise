#include <pybind11/pybind11.h>
#include "../../include/compartments/Synapse.hpp"

namespace py = pybind11;

void bind_Synapse(py::module& m) {
  py::class_<Synapse, Compartment>(m, "Synapse")
    .def(py::init<Compartment&, const std::string&>());
}
