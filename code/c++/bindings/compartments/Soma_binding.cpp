#include <pybind11/pybind11.h>
#include "../../include/compartments/Soma.hpp"

namespace py = pybind11;

void bind_Soma(py::module& m) {
  py::class_<Soma, Compartment>(m, "Soma")
    .def(py::init<const std::string&, const double&>(),
         py::arg("name")="no_name_Soma", py::arg("length")=10,
         "Create soma with certain name and length")
    .def(py::init<const std::string&, const double&, const double&, const double&, const double&, const double&, const size_t&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&>(),
         py::arg("name")="no_name_Soma", py::arg("length")=10, py::arg("x")=0, py::arg("y")=0, py::arg("z")=0, py::arg("radius")=5, py::arg("n_gene_copies")=1, py::arg("gene_activation_rate")=1/12., py::arg("gene_deactivation_rate")=1/12., py::arg("mRNA_decay_rate")=0.0432, py::arg("translation_rate")=75.6, py::arg("protein_decay_rate")=0.004356, py::arg("mRNA_diffusion_constant")=3.4e-3, py::arg("protein_diffusion_constant")=.24, py::arg("mRNA_forward_trafficking_velocity")=.5e-2, py::arg("mRNA_backward_trafficking_velocity")=.1e-2, py::arg("protein_forward_trafficking_velocity")=0, py::arg("protein_backward_trafficking_velocity")=0,
         "Create soma with certain name and length, location (x,y,z) and radius");
}
