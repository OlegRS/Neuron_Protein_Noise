#ifndef __ANALYTIC_ENGINE_HPP__
#define __ANALYTIC_ENGINE_HPP__

#include <armadillo>

#include "../../Neuron.hpp"
#include "../../compartments/Soma.hpp"
#include "../../compartments/Dendritic_segment.hpp"
#include "../../compartments/Synapse.hpp"
#include "../../junctions/Som_den_junction.hpp"
#include "../../junctions/Den_syn_junction.hpp"
#include "../../junctions/Den_den_junction.hpp"

class Analytic_engine {
  // Parameters
  Neuron* p_neuron = NULL;
  unsigned int o1_dim, o2_dim;
  arma::mat *p_o1_mat, *p_o2_mat,
    o1_mRNA_matrix, o1_prot_matrix,
    *o2_gene_mRNA_mat, *o2_gene_prot_mat,
    *o2_mRNA_mRNA_mat, *o2_mRNA_prot_mat, *o2_prot_prot_mat,
    *o2_nonstationary_RHS_mat;
  arma::vec *p_o1_RHS, *p_o2_RHS,
    expectations, *p_covariances,
    o1_mRNA_RHS, o1_prot_RHS,
    mRNA_expectations, protein_expectations,
    *o2_gene_mRNA, *o2_gene_mRNA_RHS,
    *o2_gene_prot, *o2_gene_prot_RHS,
    *o2_mRNA_mRNA, *o2_mRNA_mRNA_RHS,
    *o2_mRNA_prot, *o2_mRNA_prot_RHS,
    *o2_prot_prot, *o2_prot_prot_RHS;

  // Variables
  std::vector<std::string> o1_var_names, o1_mRNA_names, o1_prot_names,  *p_o2_var_names;
  std::vector<double*> p_o1_vars,
    p_o1_mRNA_expectations, p_o1_prot_expectations; // Pointers to variables within the compartments
  // /// These are needed to fill o1_matrix
  // unsigned int parent_start_ind = 0;
  // unsigned int desc_start_ind = 3;

  Analytic_engine& initialise_o1_mat_and_RHS();
  // Sets 1st order matrix starting from the given compartment
  const Compartment* set_o1_soma();
  void set_o1_matrix(const Compartment&);
  const Compartment* sem_set_o1_soma();
  void sem_set_o1_matrix(const Compartment&);
  Analytic_engine& internalise_expectations(); // Writes expectations into compartments

  // 1st order with separations of gene-mRNA-protein dynamics
  const Compartment* set_o1_mRNA_soma(); //Computes expectation of active genes and sets RHS for mRNA eqns
  const Compartment* set_o1_prot_soma();
  void set_o1_mRNA_matrix(const Compartment&);
  Analytic_engine& internalise_mRNA_expectations();
  Analytic_engine& internalise_prot_expectations();
  void set_o1_prot_matrix(const Compartment&);

  // 2nd order computation
  void initialise_o2();
  void set_prot_index_from(Compartment& compartment); // Needed for sem_o2
  void sem_initialise_o2();
  const Compartment* sem_set_soma();
  void sem_set_expectations(const Compartment& parent);
  void set_o2_soma();
  void set_o2_matrix();
  void set_o2_nonstationary_RHS_soma();
  void sem_set_o2_nonstationary_RHS_soma();
  void set_o2_nonstationary_RHS_mat();
  void sem_set_o2_nonstationary_RHS_mat();

  void sem_set_o2_soma();
  void sem_set_o2_matrix();

  double o2_gene_gene;
  const Compartment* set_o2_gene_mRNA_soma();
  void set_o2_gene_mRNA_matrix(const Compartment&);
  const Compartment* set_o2_gene_prot_soma();
  void set_o2_gene_prot_matrix(const Compartment&);
  const Compartment* set_o2_mRNA_mRNA_soma();
  void set_o2_mRNA_mRNA_matrix(const Compartment&);
  const Compartment* set_o2_mRNA_prot_soma();
  void set_o2_mRNA_prot_matrix(const Compartment&);
  const Compartment* set_o2_prot_prot_soma();
  void set_o2_prot_prot_matrix(const Compartment&);
   
  void set_o2_RHS();
  void sem_set_o2_RHS();
  Analytic_engine& clear_o1_mat_and_RHS();
  Analytic_engine& clear_o1();
  Analytic_engine& clear_o2_mat_and_RHS();
  Analytic_engine& clear_o2();
    
public:

  size_t o2_ind(const size_t &i, const size_t &j, const size_t &dim) const; // o2 index conversion
  size_t o2_ind(const size_t &i, const size_t &j) const {return o2_ind(i,j,o1_dim);}
  size_t o2_ind_asym(const size_t &i, const size_t &j, const size_t &dim_x) const {return dim_x*i+j;}
  inline size_t sem_o2_ind(const size_t &i, const size_t &j) const; // semantic o2 index conversion

  
  Analytic_engine(Neuron& neuron) :
    p_neuron(&neuron),
    o1_dim(3 + 2*neuron.p_dend_segments.size() + neuron.p_synapses.size()),
    o2_dim(o1_dim*(o1_dim+1)/2),
    p_o1_mat(NULL), p_o1_RHS(NULL),
    o1_mRNA_matrix(1+neuron.p_dend_segments.size(), 1+neuron.p_dend_segments.size()),
    o1_prot_matrix(1+neuron.p_dend_segments.size()+neuron.p_synapses.size(), 1+neuron.p_dend_segments.size()+neuron.p_synapses.size()),
    o1_mRNA_RHS(1+neuron.p_dend_segments.size()),
    mRNA_expectations(1+neuron.p_dend_segments.size()),
    o1_mRNA_names(1+neuron.p_dend_segments.size()),
    o1_prot_RHS(1+neuron.p_dend_segments.size()+neuron.p_synapses.size()),
    protein_expectations(1+neuron.p_dend_segments.size()+neuron.p_synapses.size()),
    o1_prot_names(1+neuron.p_dend_segments.size()+neuron.p_synapses.size()),
    p_o2_mat(NULL), p_o2_RHS(NULL), p_o2_var_names(NULL), p_covariances(NULL),
    expectations(o1_dim),
    o1_var_names(o1_dim),
    p_o1_vars(o1_dim),
    p_o1_mRNA_expectations(1+neuron.p_dend_segments.size()),
    p_o1_prot_expectations(1+neuron.p_dend_segments.size()+neuron.p_synapses.size()),
    o2_gene_gene(neuron.p_soma->gene_activation_rate*(neuron.p_soma->number_of_gene_copies-1)/(neuron.p_soma->gene_activation_rate + neuron.p_soma->gene_deactivation_rate)*neuron.p_soma->n_active_genes_expectation),
    o2_gene_mRNA(NULL), o2_gene_mRNA_RHS(NULL), o2_gene_mRNA_mat(NULL),
    o2_gene_prot(NULL), o2_gene_prot_RHS(NULL), o2_gene_prot_mat(NULL),
    o2_mRNA_mRNA(NULL), o2_mRNA_mRNA_RHS(NULL), o2_mRNA_mRNA_mat(NULL),
    o2_mRNA_prot(NULL), o2_mRNA_prot_RHS(NULL), o2_mRNA_prot_mat(NULL),
    o2_prot_prot(NULL), o2_prot_prot_RHS(NULL), o2_prot_prot_mat(NULL),
    o2_nonstationary_RHS_mat(NULL)
  {neuron.prot_ind = neuron.p_dend_segments.size()+1;}
  
  std::string master_equation() const; // Not implemented
  std::string generating_function_PDE() const; // Not implemented

  Analytic_engine& initialize();

  Analytic_engine& stationary_expectations();//const Neuron& neur = *p_neuron);
  Analytic_engine& sem_stationary_expectations();//const Neuron& neur = *p_neuron);
  Analytic_engine& mRNA_stationary_expectations();
  Analytic_engine& mRNA_o1_eigen_decomposition();
  Analytic_engine& protein_stationary_expectations();
  Analytic_engine& protein_o1_eigen_decomposition();
  Analytic_engine& nonstationary_expectations(const std::list<double>& times);//const Neuron& neur = *p_neuron);
  Analytic_engine& sem_nonstationary_expectations(const std::list<double>& times);
  
  Analytic_engine& stationary_covariances();//const Neuron& neur = *p_neuron);
  Analytic_engine& sem_stationary_covariances(); // Semantically grouped stationary covariances
  Analytic_engine& sem_stationary_pearson_correlations();

  Analytic_engine& gene_mRNA_stationary_covariances();
  Analytic_engine& gene_protein_stationary_covariances();
  Analytic_engine& mRNA_mRNA_stationary_covariances();
  Analytic_engine& mRNA_protein_stationary_covariances();
  Analytic_engine& protein_protein_stationary_covariances();
  Analytic_engine& nonstationary_covariances(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2);
  Analytic_engine& sem_nonstationary_covariances(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2);
  Analytic_engine& sem_nonstationary_covariances_using_integral(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2);
  Analytic_engine& sem_nonstationary_covariances_direct_ODE_solver(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2);

  arma::vec* G1() {return &expectations;}
  std::vector<std::string>* o1_variables_names() {return &o1_var_names;}
  arma::vec* G2() {return p_covariances;}

  size_t dim_o1() const {return o1_dim;}
  
  // std::vector<std::pair<std::string,double>> expectations();
  // std::vector<std::vector<std::pair<std::string,double>>> Analytic_engine& correlations();

  Analytic_engine& clear_all();
    
  ~Analytic_engine();
};

#endif
