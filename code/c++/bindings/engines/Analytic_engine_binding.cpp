#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../include/engines/analytic/Analytic_engine.hpp"

namespace py = pybind11;

void bind_Analytic_engine(py::module& m) {
    py::class_<Analytic_engine>(m, "Analytic_engine")
      .def(py::init<Neuron&, bool>(),
           py::arg("neuron"), py::arg("cov_mat_init")=false,
           "Initialise Analytic Engine for a given neuron")
      .def("stationary_expectations_and_correlations", &Analytic_engine::stationary_expectations_and_correlations,
           "Compute all stationary expectations and correlations")
      .def("stationary_mRNA_expectations", &Analytic_engine::stationary_mRNA_expectations,
           "Compute stationary expectations of mRNA counts using matrix inversion")
      .def("stationary_protein_expectations", &Analytic_engine::stationary_protein_expectations,
           "Compute stationary expectations of protein counts using matrix inversion.\n!!! THIS CAN ONLY BE CALLED WHEN mRNA EXPECTATIONS ARE COMPUTED, SO CALL stationary_mRNA_expectations() FIRST !!!")
      .def("stationary_gene_mRNA_covariances", &Analytic_engine::stationary_gene_mRNA_covariances,
           "Compute correlations (<nm>) of active gene counts with protein counts using matrix inversion.\n!!! THIS CAN ONLY BE CALLED WHEN mRNA EXPECTATIONS ARE COMPUTED, SO CALL stationary_mRNA_expectations() FIRST !!!")
      .def("stationary_gene_protein_covariances", &Analytic_engine::stationary_gene_protein_covariances,
           "Compute stationary correlations (<n*m>) of active gene counts with protein counts using matrix inversion.\n!!! THIS CAN ONLY BE CALLED WHEN mRNA EXPECTATIONS and gene-mRNA CORRELATIONS ARE COMPUTED, SO CALL stationary_mRNA_expectations() and stationary_gene_mRNA_covariances() FIRST !!!")
      .def("stationary_mRNA_mRNA_covariances", &Analytic_engine::stationary_mRNA_mRNA_covariances,
           "Compute stationary correlations (<m*m>) of active mRNA counts with mRNA counts using matrix inversion.\n!!! THIS CAN ONLY BE CALLED WHEN mRNA EXPECTATIONS and gene-mRNA CORRELATIONS ARE COMPUTED, SO CALL stationary_mRNA_expectations() and stationary_gene_mRNA_covariances() FIRST !!!")
      .def("stationary_mRNA_protein_covariances", &Analytic_engine::stationary_mRNA_protein_covariances,
           "Compute stationary correlations (<m*p>) of active mRNA counts with protein counts using matrix inversion.\n!!! THIS CAN ONLY BE CALLED WHEN mRNA EXPECTATIONS and gene-mRNA CORRELATIONS ARE COMPUTED, SO CALL stationary_mRNA_expectations() and stationary_gene_mRNA_covariances() FIRST !!!")
      .def("stationary_protein_protein_covariances", &Analytic_engine::stationary_protein_protein_covariances,
           "Compute stationary correlations (<p*p>) of active protein counts with protein counts using matrix inversion.\n!!! THIS CAN ONLY BE CALLED WHEN mRNA EXPECTATIONS and gene-mRNA CORRELATIONS ARE COMPUTED, SO CALL stationary_mRNA_expectations() and stationary_gene_mRNA_covariances() FIRST !!!")
      .def("active_genes_expectation", &Analytic_engine::active_genes_expectation,
           "Return expected number of active gene copies")
      .def("mRNA_expectation", &Analytic_engine::mRNA_expectation,
           py::arg("compartment"),
           "Return expected number of mRNAs in the compartment")
      .def("protein_expectation", &Analytic_engine::protein_expectation,
           py::arg("compartment"),
           "Return expected number of proteins in the compartment")
      .def("gene_mRNA_correlation", &Analytic_engine::gene_mRNA_correlation,
           py::arg("compartment"),
           "Return gene-mRNA correlation <n_g*n_mRNA>")
      .def("mRNA_mRNA_correlation", &Analytic_engine::mRNA_mRNA_correlation,
           py::arg("compartment"), py::arg("compartment"),
           "Return mRNA-mRNA correlation <n_mRNA*n_mRNA>")
      .def("gene_protein_correlation", &Analytic_engine::gene_protein_correlation,
           py::arg("compartment"),
           "Return gene-protein correlation <n_g*n_mRNA>")
      .def("mRNA_protein_correlation", &Analytic_engine::mRNA_protein_correlation,
           py::arg("compartment"), py::arg("compartment"),
           "Return mRNA-protein correlation <n_mRNA*n_prot>")
      .def("protein_protein_correlation", &Analytic_engine::protein_protein_correlation,
           py::arg("compartment"), py::arg("compartment"),
           "Return protein-protein correlation <n_prot*n_prot>");
}
