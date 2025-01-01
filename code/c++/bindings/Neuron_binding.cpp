#include <pybind11/pybind11.h>
#include <sstream>
#include "../include/Neuron.hpp"
#include "../include/compartments/Soma.hpp"

namespace py = pybind11;

// Helper function to invoke std::ostream <<
std::string to_string(const Neuron& neuron) {
    std::ostringstream oss;
    oss << neuron;  // Uses the overloaded << operator
    return oss.str();
}

void bind_Neuron(py::module& m) {    
  py::class_<Neuron>(m, "Neuron")
    .def(py::init<Soma&, const std::string&>())
    .def("__str__", &to_string); // Link Python's str() to C++'s << 
}
