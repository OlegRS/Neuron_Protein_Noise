#include <pybind11/pybind11.h>
#include "../include/Neuron.hpp"

namespace py = pybind11;

void bind_Neuron(py::module& m) {
  // // Declare Soma as an opaque type
  // py::class_<Soma>(m, "Soma", py::module_local());
    
  py::class_<Neuron>(m, "Neuron")
    .def(py::init<Soma&, const std::string&>());
}
