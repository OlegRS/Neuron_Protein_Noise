#include <pybind11/pybind11.h>
#include "../../include/engines/analytic/Analytic_engine.hpp"

namespace py = pybind11;

void bind_Analytic_engine(py::module& m) {
    py::class_<Analytic_engine>(m, "Analytic_engine")
      .def(py::init<Neuron&, bool>())
      .def("stationary_expectations_and_correlations", &Analytic_engine::stationary_expectations_and_correlations);
}
