#include "Analytic_engine.hpp"
#include <math.h>

const Compartment* Analytic_engine::set_o1_mRNA_soma() {
  auto& soma = *p_neuron->p_soma;
  o1_mRNA_matrix(0,0) -= soma.mRNA_decay_rate;
  o1_mRNA_RHS(0) = -soma.n_active_genes_expectation * soma.transcription_rate;

  o1_mRNA_names[0] = soma.name;
  p_o1_mRNA_expectations[0] = &soma.n_mRNA_expectation;

  return p_neuron->p_soma;
}

void Analytic_engine::set_o1_mRNA_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty()) {
    return;
  }

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t parent_start_ind = p_junc -> p_from -> mRNA_ind-1;
    size_t desc_start_ind = p_junc -> p_to -> mRNA_ind-1;
    
    if(p_junc->type() != DEN_SYN) {

      o1_mRNA_names[desc_start_ind] = p_junc->p_to->name;
      
      p_o1_mRNA_expectations[desc_start_ind] = &p_junc->p_to->n_mRNA_expectation;

      o1_mRNA_matrix(desc_start_ind, desc_start_ind) -= p_junc->p_to->mRNA_decay_rate;

      o1_mRNA_matrix(parent_start_ind, parent_start_ind) -= p_junc->fwd_mRNA_hop_rate;
      o1_mRNA_matrix(parent_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      o1_mRNA_matrix(desc_start_ind, parent_start_ind) += p_junc->fwd_mRNA_hop_rate;
      o1_mRNA_matrix(desc_start_ind, desc_start_ind) -= p_junc->bkwd_mRNA_hop_rate;
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o1_mRNA_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::mRNA_stationary_expectations() {

  std::cerr << "Setting the mRNA matrix...\n";
  set_o1_mRNA_matrix(*set_o1_mRNA_soma());
  std::cerr << "Inverting the mRNA matrix...\n";
  mRNA_expectations = o1_mRNA_matrix.i()*o1_mRNA_RHS;

  std::cout << "n_active_genes_expectation = " << p_neuron->p_soma->n_active_genes_expectation << '\n';
  std::cout << "mRNA_expectations:\n";
  for(size_t i=0; i<1+p_neuron->p_dend_segments.size(); ++i)
    std::cerr << o1_mRNA_names[i] << ": " << mRNA_expectations(i) << std::endl;
  
  return internalise_mRNA_expectations();
}

Analytic_engine& Analytic_engine::mRNA_o1_eigen_decomposition() {
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(1+p_neuron->p_dend_segments.size());
  arma::cx_mat eigvec_c;
  arma::mat trans(1+p_neuron->p_dend_segments.size(), 1+p_neuron->p_dend_segments.size()); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mRNA_matrix);

  std::cout << "mRNA_eigval:\n";
  for(size_t i=0; i<1+p_neuron->p_dend_segments.size(); ++i) {
    std::cout << (eigval(i) = -eigval_c(i).real()) << ',';
    for(size_t j=0; j<1+p_neuron->p_dend_segments.size(); ++j)
      trans(i,j) = -eigvec_c(i,j).real();
  }
  std::cout << "\nmRNA_eigvec:\n" << trans << std::endl;

  std::cout << "Compartment names for mRNA:\n";
  for(auto& name : o1_mRNA_names)
    std::cout << name + "__mRNA" << std::endl;

  return *this;
}

Analytic_engine& Analytic_engine::internalise_mRNA_expectations() {
  size_t mRNA_size = 1 + p_neuron->p_dend_segments.size();
  for(size_t i=0; i<mRNA_size; ++i)
    *p_o1_mRNA_expectations[i] = mRNA_expectations(i);
  return *this;
}

const Compartment* Analytic_engine::set_o1_prot_soma() {
  auto& soma = *p_neuron->p_soma;
  o1_prot_matrix(0,0) -= soma.protein_decay_rate;
  o1_prot_RHS(0) = -mRNA_expectations(0)*soma.translation_rate;

  o1_prot_names[0] = soma.name;
  p_o1_prot_expectations[0] = &soma.n_prot_expectation;

  return p_neuron->p_soma;
}

void Analytic_engine::set_o1_prot_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty()) {
    return;
  }

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t parent_start_ind = p_junc -> p_from -> id;
    size_t desc_start_ind = p_junc -> p_to -> id;
    
    o1_prot_names[desc_start_ind] = p_junc->p_to->name;      
    p_o1_prot_expectations[desc_start_ind] = &p_junc->p_to->n_prot_expectation;

    o1_prot_matrix(parent_start_ind, parent_start_ind) -= p_junc->fwd_prot_hop_rate;
    o1_prot_matrix(parent_start_ind, desc_start_ind) += p_junc->bkwd_prot_hop_rate;
    o1_prot_matrix(desc_start_ind, parent_start_ind) += p_junc->fwd_prot_hop_rate;
    o1_prot_matrix(desc_start_ind, desc_start_ind) -= p_junc->bkwd_prot_hop_rate;


    if(p_junc->type() != DEN_SYN) {
      o1_prot_matrix(desc_start_ind, desc_start_ind) -= p_junc->p_to->protein_decay_rate;
      o1_prot_RHS(desc_start_ind) = -(p_junc->p_to->n_mRNA_expectation)*(p_junc->p_to->translation_rate);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o1_prot_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::protein_stationary_expectations() {

  std::cerr << "Setting the prot matrix...\n";
  set_o1_prot_matrix(*set_o1_prot_soma());
  std::cerr << "Inverting the prot matrix...\n";
  protein_expectations = o1_prot_matrix.i()*o1_prot_RHS;

  std::cout << "Protein_expectations:\n";
  for(size_t i=0; i<1+p_neuron->p_dend_segments.size()+p_neuron->p_synapses.size(); ++i)
    std::cerr << o1_prot_names[i] << ": " << protein_expectations(i) << std::endl;
  
  return internalise_prot_expectations();
}

Analytic_engine& Analytic_engine::protein_o1_eigen_decomposition() {
  std::cerr << "Computing eigen decomposition...\n";
  size_t sz = 1+p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size();
  arma::cx_vec eigval_c;
  arma::vec eigval(sz);
  arma::cx_mat eigvec_c;
  arma::mat trans(sz, sz); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_prot_matrix);

  std::cout << "Protein_eigval:\n";
  for(size_t i=0; i<sz; ++i) {
    std::cout << (eigval(i) = -eigval_c(i).real()) << ',';
    for(size_t j=0; j<sz; ++j)
      trans(i,j) = -eigvec_c(i,j).real();
  }
  std::cout << "\nProtein_eigvec:\n" << trans << std::endl;

  std::cout << "Compartment names for proteins:\n";
  for(auto& name : o1_prot_names)
    std::cout << name + "__Prot" << std::endl;

  return *this;
}

Analytic_engine& Analytic_engine::internalise_prot_expectations() {
  size_t prot_size = 1 + p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size();
  for(size_t i=0; i<prot_size; ++i)
    *p_o1_prot_expectations[i] = protein_expectations(i);
  return *this;
}

const Compartment* Analytic_engine::set_o1_soma() {
  
  auto& soma = *p_neuron->p_soma;

  o1_mat(0,0) = -soma.gene_activation_rate - soma.gene_deactivation_rate;

  o1_mat(1,0) = soma.transcription_rate;
  o1_mat(1,1) = -soma.mRNA_decay_rate;

  o1_mat(2,1) = soma.translation_rate;
  o1_mat(2,2) = -soma.protein_decay_rate;

  o1_var_names[0] = soma.name + "__Gene";
  o1_var_names[1] = soma.name + "__mRNA";
  o1_var_names[2] = soma.name + "__Prot";

  p_o1_vars[0] = &soma.n_active_genes_expectation;
  p_o1_vars[1] = &soma.n_mRNA_expectation;
  p_o1_vars[2] = &soma.n_prot_expectation;

  // return static_cast<const Compartment*&>(neur.p_soma);
  return p_neuron->p_soma;
}

void Analytic_engine::set_o1_matrix(const Compartment& parent) {
  
  if(parent.it_p_out_junctions.empty()) {
    return;
  }

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t& parent_start_ind = p_junc -> p_from -> o1_index;
    size_t& desc_start_ind = p_junc -> p_to -> o1_index;
    
    if(p_junc->type() == DEN_SYN) {
      // o1_var_names.push_back(p_junc->p_to->name + "__Prot");
      o1_var_names[p_junc->p_to->o1_index] = p_junc->p_to->name + "__Prot";
      p_o1_vars[p_junc->p_to->o1_index] = &p_junc->p_to->n_prot_expectation;
      
      o1_mat(parent_start_ind+1, parent_start_ind+1) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_start_ind+1, desc_start_ind) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_start_ind, parent_start_ind+1) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == DEN_DEN) {
      // o1_var_names.push_back(p_junc->p_to->name + "__mRNA");
      // o1_var_names.push_back(p_junc->p_to->name + "__Prot");

      o1_var_names[p_junc->p_to->o1_index] = p_junc->p_to->name + "__mRNA";
      o1_var_names[p_junc->p_to->o1_index+1] = p_junc->p_to->name + "__Prot";
      p_o1_vars[p_junc->p_to->o1_index] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[p_junc->p_to->o1_index+1] = &p_junc->p_to->n_prot_expectation;

      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->p_to->mRNA_decay_rate;
      o1_mat(desc_start_ind+1, desc_start_ind) += p_junc->p_to->translation_rate;
      o1_mat(desc_start_ind+1, desc_start_ind+1) -= p_junc->p_to->protein_decay_rate;

      o1_mat(parent_start_ind, parent_start_ind) -= p_junc->fwd_mRNA_hop_rate;
      o1_mat(parent_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      o1_mat(desc_start_ind, parent_start_ind) += p_junc->fwd_mRNA_hop_rate;
      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->bkwd_mRNA_hop_rate;

      o1_mat(parent_start_ind+1, parent_start_ind+1) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_start_ind+1, desc_start_ind+1) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_start_ind+1, parent_start_ind+1) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_start_ind+1, desc_start_ind+1) -= p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == SOM_DEN) {
      o1_var_names[p_junc->p_to->o1_index] = p_junc->p_to->name + "__mRNA";
      o1_var_names[p_junc->p_to->o1_index+1] = p_junc->p_to->name + "__Prot";
      p_o1_vars[p_junc->p_to->o1_index] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[p_junc->p_to->o1_index+1] = &p_junc->p_to->n_prot_expectation;

      // o1_var_names.push_back(p_junc->p_to->name + "__mRNA");
      // o1_var_names.push_back(p_junc->p_to->name + "__Prot");
      
      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->p_to->mRNA_decay_rate;
      o1_mat(desc_start_ind+1, desc_start_ind) += p_junc->p_to->translation_rate;
      o1_mat(desc_start_ind+1, desc_start_ind+1) -= p_junc->p_to->protein_decay_rate;

      o1_mat(parent_start_ind+1, parent_start_ind+1) -= p_junc->fwd_mRNA_hop_rate;
      o1_mat(parent_start_ind+1, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      o1_mat(desc_start_ind, parent_start_ind+1) += p_junc->fwd_mRNA_hop_rate;
      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->bkwd_mRNA_hop_rate;

      o1_mat(parent_start_ind+2, parent_start_ind+2) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_start_ind+2, desc_start_ind+1) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_start_ind+1, parent_start_ind+2) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_start_ind+1, desc_start_ind+1) -= p_junc->bkwd_prot_hop_rate;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o1_matrix(*((*it_p_junc)->p_to));
}

const Compartment* Analytic_engine::sem_set_o1_soma() {
  
  auto& soma = *p_neuron->p_soma;

  set_prot_index_from(soma);

  o1_mat(0,0) = -soma.gene_activation_rate - soma.gene_deactivation_rate;

  o1_mat(1,0) = soma.transcription_rate; //soma.mRNA_ind==1
  o1_mat(1,1) = -soma.mRNA_decay_rate;

  o1_mat(soma.prot_ind,1) = soma.translation_rate;
  o1_mat(soma.prot_ind,soma.prot_ind) = -soma.protein_decay_rate;

  o1_var_names[0] = soma.name + "__Gene";
  o1_var_names[1] = soma.name + "__mRNA";
  o1_var_names[soma.prot_ind] = soma.name + "__Prot";

  p_o1_vars[0] = &soma.n_active_genes_expectation;
  p_o1_vars[1] = &soma.n_mRNA_expectation;
  p_o1_vars[soma.prot_ind] = &soma.n_prot_expectation;

  // return static_cast<const Compartment*&>(neur.p_soma);
  return p_neuron->p_soma;
}

void Analytic_engine::sem_set_o1_matrix(const Compartment& parent) {

  if(parent.it_p_out_junctions.empty()) {
    return;
  }

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t& parent_prot_ind = p_junc -> p_from -> prot_ind;
    size_t& desc_prot_ind = p_junc -> p_to -> prot_ind;
    size_t& parent_mRNA_ind = p_junc -> p_from -> mRNA_ind;
    size_t& desc_mRNA_ind = p_junc -> p_to -> mRNA_ind;

    if(p_junc->type() == DEN_SYN) {
      // o1_var_names.push_back(p_junc->p_to->name + "__Prot");
      o1_var_names[desc_prot_ind] = p_junc->p_to->name + "__Prot";
      p_o1_vars[desc_prot_ind] = &p_junc->p_to->n_prot_expectation;
      
      o1_mat(parent_prot_ind, parent_prot_ind) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == DEN_DEN) {
      // o1_var_names.push_back(p_junc->p_to->name + "__mRNA");
      // o1_var_names.push_back(p_junc->p_to->name + "__Prot");

      o1_var_names[desc_mRNA_ind] = p_junc->p_to->name + "__mRNA";
      o1_var_names[desc_prot_ind] = p_junc->p_to->name + "__Prot";
      p_o1_vars[desc_mRNA_ind] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[desc_prot_ind] = &p_junc->p_to->n_prot_expectation;

      o1_mat(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->p_to->mRNA_decay_rate;
      o1_mat(desc_prot_ind, desc_mRNA_ind) += p_junc->p_to->translation_rate;
      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->p_to->protein_decay_rate;

      o1_mat(parent_mRNA_ind, parent_mRNA_ind) -= p_junc->fwd_mRNA_hop_rate;
      o1_mat(parent_mRNA_ind, desc_mRNA_ind) += p_junc->bkwd_mRNA_hop_rate;
      o1_mat(desc_mRNA_ind, parent_mRNA_ind) += p_junc->fwd_mRNA_hop_rate;
      o1_mat(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->bkwd_mRNA_hop_rate;

      o1_mat(parent_prot_ind, parent_prot_ind) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == SOM_DEN) {
      o1_var_names[desc_mRNA_ind] = p_junc->p_to->name + "__mRNA";
      o1_var_names[desc_prot_ind] = p_junc->p_to->name + "__Prot";
      p_o1_vars[desc_mRNA_ind] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[desc_prot_ind] = &p_junc->p_to->n_prot_expectation;
      
      o1_mat(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->p_to->mRNA_decay_rate;
      o1_mat(desc_prot_ind, desc_mRNA_ind) += p_junc->p_to->translation_rate;
      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->p_to->protein_decay_rate;

      o1_mat(parent_mRNA_ind, parent_mRNA_ind) -= p_junc->fwd_mRNA_hop_rate;
      o1_mat(parent_mRNA_ind, desc_mRNA_ind) += p_junc->bkwd_mRNA_hop_rate;
      o1_mat(desc_mRNA_ind, parent_mRNA_ind) += p_junc->fwd_mRNA_hop_rate;
      o1_mat(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->bkwd_mRNA_hop_rate;

      o1_mat(parent_prot_ind, parent_prot_ind) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->bkwd_prot_hop_rate;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    sem_set_o1_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::internalise_expectations() {
  for(size_t i=0; i<o1_dim; ++i)
    *p_o1_vars[i] = expectations(i);
  return *this;
}

Analytic_engine& Analytic_engine::stationary_expectations() {
  o1_RHS(0) = -(p_neuron->p_soma->gene_activation_rate) * (p_neuron->p_soma->number_of_gene_copies);
  std::cerr << "Setting o1_matrix...\n";
  set_o1_matrix(*set_o1_soma());
  std::cerr << "Inverting o1_matrix...\n";
  arma::mat inv_o1_matrix = o1_mat.i();
  std::cerr << "Done with o1_matrix inversion\n";

  expectations = inv_o1_matrix * o1_RHS;

  size_t i=0;
  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << "=\t" << expectations[i++] << std::endl;
  
  return internalise_expectations();
}

Analytic_engine& Analytic_engine::sem_stationary_expectations() {
  o1_RHS(0) = -(p_neuron->p_soma->gene_activation_rate) * (p_neuron->p_soma->number_of_gene_copies);
  std::cerr << "Setting o1_matrix...\n";
  sem_set_o1_matrix(*sem_set_o1_soma());
  std::cerr << "Inverting o1_matrix...\n";
  arma::mat inv_o1_matrix = o1_mat.i();
  std::cerr << "Done with o1_matrix inversion\n";

  expectations = inv_o1_matrix * o1_RHS;

  size_t i=0;
  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << "=\t" << expectations[i++] << std::endl;
  
  return internalise_expectations();
}

Analytic_engine& Analytic_engine::nonstationary_expectations(const std::list<double>& times) { //const Neuron& neur) {
  o1_RHS(0) = -(p_neuron->p_soma->gene_activation_rate) * (p_neuron->p_soma->number_of_gene_copies);
  std::cerr << "Setting o1_matrix...\n";
  set_o1_matrix(*set_o1_soma());
  o1_mat = -o1_mat;
  
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;

  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    std::cout << "stationary_expectations[" << i << "]=" << (stationary_expectations(i) = (p_neuron->p_soma->gene_activation_rate)*(p_neuron->p_soma->number_of_gene_copies)*sum)  << std::endl;
  }

  // Setting integrals of motion c = G(0) - A*b, where G(0)=0 for now
  arma::vec c = 0 - stationary_expectations;

  for(auto& t : times) {
    std::cout << t << ',';
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k) {
        double s_sum=0;
        for(size_t s=0; s<o1_dim; ++s)
          s_sum += inv_tm(k,s)*c(s);
        k_sum += exp(-eigval(k)*t) * tm(i,k) * s_sum;
      }
      expectations(i) = stationary_expectations(i) + k_sum;
      std::cout << expectations(i) << ',';
    }
    std::cout << std::endl;
  }

  return internalise_expectations();
}


Analytic_engine& Analytic_engine::sem_nonstationary_expectations(const std::list<double>& times) { //const Neuron& neur) {
  o1_RHS(0) = -(p_neuron->p_soma->gene_activation_rate) * (p_neuron->p_soma->number_of_gene_copies);
  std::cerr << "Setting o1_matrix...\n";
  o1_mat.zeros();
  sem_set_o1_matrix(*sem_set_o1_soma());
  o1_mat = -o1_mat;

  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;
    
  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    std::cout << "stationary_expectations[" << i << "]=" << (stationary_expectations(i) = (p_neuron->p_soma->gene_activation_rate)*(p_neuron->p_soma->number_of_gene_copies)*sum)  << std::endl;
  }

  // Setting integrals of motion c = G(0) - A*b, where G(0)=0 for now
  arma::vec c = 0 - stationary_expectations;

  for(auto& t : times) {
    std::cout << t << ',';
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k) {
        double s_sum=0;
        for(size_t s=0; s<o1_dim; ++s)
          s_sum += inv_tm(k,s)*c(s);
        k_sum += exp(-eigval(k)*t) * tm(i,k) * s_sum;
      }
      expectations(i) = stationary_expectations(i) + k_sum;
      std::cout << expectations(i) << ',';
    }
    std::cout << std::endl;
  }

  return internalise_expectations();
}

Analytic_engine& Analytic_engine::nonstationary_covariances(const std::list<double>& times) {
  o1_RHS(0) = -(p_neuron->p_soma->gene_activation_rate) * (p_neuron->p_soma->number_of_gene_copies);
  std::cerr << "Setting o1_matrix...\n";
  o1_mat.zeros();
  set_o1_matrix(*set_o1_soma());
  o1_mat = -o1_mat;
  
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;

  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    stationary_expectations(i) = p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * sum;
  }

  for(size_t i=0; i<o1_dim; ++i)
    std::cout << o1_var_names[i] << ": " << stationary_expectations(i) << std::endl;

  // Setting integrals of motion c = G(0) - A*b, where G(0)=0 for now
  arma::vec c_vec = 0 - stationary_expectations;

  // Initialising second order
  initialise_o2();
  set_o2_matrix();
  set_o2_nonstationary_RHS_mat();

  // Debug_start
  std::cerr << "Nonstationary RHS_mat:\n" << *o2_nonstationary_RHS_mat << std::endl;
  // this->stationary_expectations();
  set_o2_RHS();
  std::cerr << "(*p_o2_RHS):\n" << (*p_o2_RHS) << std::endl;
  for(size_t j=0; j<o2_dim; ++j) {
    double sum = 0;
    for(size_t i=0; i<o1_dim; ++i)
      sum -= (*o2_nonstationary_RHS_mat)(j,i)*stationary_expectations(i);
    std::cerr << "sum = " << sum << ", (*p_o2_RHS)(j) = " << (*p_o2_RHS)(j) << std::endl;
  }
  std::cerr << std::endl;
  // Debug_end

  std::vector<double> k_sum(o1_dim); // to store precomputed k-sum
  
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  o2_mat = -o2_mat;
  
  std::cerr << "Computing o2 eigen decomposition...\n";
  arma::cx_vec o2_eigval_c;
  arma::vec o2_eigval(o2_dim);
  arma::cx_mat o2_eigvec_c;
  arma::mat o2_tm(o2_dim, o2_dim); // Transition matrix
  arma::eig_gen(o2_eigval_c, o2_eigvec_c, o2_mat);

  for(size_t i=0; i<o2_dim; ++i) {
    o2_eigval(i) = o2_eigval_c(i).real();
    for(size_t j=0; j<o2_dim; ++j)
      o2_tm(i,j) = o2_eigvec_c(i,j).real();
  }

  std::cerr << "o2_eigenvalues:\n"
            << o2_eigval << std::endl;
            // << "o2_eigenvectors:\n"
            // << o2_tm << std::endl;

  std::cerr << "Inverting o2 transition matrix...\n";
  arma::mat o2_inv_tm = o2_tm.i(); // Transition matrix
  std::cerr << "o2 transition matrix inverted\n";
  
  for(auto& o2_var_name : *p_o2_var_names)
    std::cout << o2_var_name << ',';
  std::cout << std::endl;

  // main loop
  std::cout << "Main loop...\n";
  for(auto& t : times) {
    std::cout << "t=" << t << '\n';
    // o1
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k) {
        double l_sum=0;
        for(size_t l=0; l<o1_dim; ++l)
          l_sum += inv_tm(k,l)*c_vec(l);
        k_sum += exp(-eigval(k)*t) * tm(i,k) * l_sum;
      }
      expectations(i) = stationary_expectations(i) + k_sum;
    }
    // o2
    for(size_t i=0; i<o2_dim; ++i) {
      (*p_covariances)[i] = 0;
      for(size_t j=0; j<o2_dim; ++j) {
        double s_sum=0;
        for(size_t s=0; s<o2_dim; ++s)
          s_sum += o2_inv_tm(j,s) * 0; // G^2_s(t=0) := 0 for now

        for(size_t zeta=0; zeta<o1_dim; ++zeta) { // Precomputing k-sum for all zeta (eta)
          k_sum[zeta] = 0;
          for(size_t k=0; k<o2_dim; ++k)
            k_sum[zeta] += o2_inv_tm(j,k)*(*o2_nonstationary_RHS_mat)(k,zeta);
        }

        double eta_sum=0;
        for(size_t eta=0; eta<o1_dim; ++eta)
          eta_sum += k_sum[eta]*stationary_expectations(eta);

        double mu_sum=0;
        for(size_t mu=0; mu<o1_dim; ++mu) {
          double zeta_sum=0;
          for(size_t zeta=0; zeta<o1_dim; ++zeta)
            zeta_sum += k_sum[zeta]*tm(zeta,mu);
          double nu_sum=0;
          for(size_t nu=0; nu<o1_dim; ++nu)
            nu_sum += inv_tm(mu,nu)*c_vec(nu);
          
          if(o2_eigval(j) != eigval(mu))
            mu_sum += (exp(-eigval(mu)*t)-exp(-o2_eigval(j)*t))/(o2_eigval(j)-eigval(mu)) * zeta_sum * nu_sum;
          else
            mu_sum += exp(-o2_eigval(j)*t) * zeta_sum * nu_sum;
            
          // else {
          //   std::cerr << "j=" << j << ", mu=" << mu << std::endl
          //             << "(*p_o2_var_names)[j] = " << (*p_o2_var_names)[j] << ", o1_var_names[mu] = " << o1_var_names[mu] << std::endl;            
          // }
        }
        (*p_covariances)[i] += o2_tm(i,j) * (exp(-o2_eigval(j)*t)*s_sum + (1-exp(-o2_eigval(j)*t))/o2_eigval(j)*eta_sum + mu_sum);
        // std::cerr << "t=" << t << std::endl;
        // std::cerr << "mu_sum[" << i << ',' << j << "] = " << (1-exp(-o2_eigval(j)*t))/o2_eigval(j)*eta_sum + mu_sum << std::endl; 
      }
      std::cerr << "t=" << t << ", " << (*p_o2_var_names)[i] + "=" << (*p_covariances)[i] << std::endl;
    }

    //////////// COMPUTING VARIANCES //////////
    std::vector<double> rmss(o1_dim);
    for(size_t i=0; i<o1_dim; ++i) {
      rmss[i] = ((*p_covariances)(o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
      if(rmss[i]>=0) {
        std::cout << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i] << std::endl;
      }
      else
        std::cout << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i] << " NEGATIVE!\n";
    }
    // std::cout << t << ", " << expectations(2) << ", " << rmss[2] << std::endl;
    //////////////////////////////////////////
  }

  return internalise_expectations();
}

Analytic_engine& Analytic_engine::sem_nonstationary_covariances(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2) {
  o1_RHS(0) = -(p_neuron->p_soma->gene_activation_rate) * (p_neuron->p_soma->number_of_gene_copies);
  std::cerr << "Setting o1_matrix...\n";
  o1_mat.zeros();
  sem_set_o1_matrix(*sem_set_o1_soma());
  o1_mat = -o1_mat;
  
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;

  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    stationary_expectations(i) = p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * sum;
  }

  std::cout << "Stationary expectations from nonstationary covariances algorithm:\n";
  for(size_t i=0; i<o1_dim; ++i)
    std::cout << o1_var_names[i] << ": " << stationary_expectations(i) << std::endl;

  // Setting integrals of motion c = G(0) - A*b
  arma::vec c_vec = *initial_G1 - stationary_expectations;

  // Initialising second order
  sem_initialise_o2();
  sem_set_o2_matrix();
  sem_set_o2_nonstationary_RHS_mat();
      
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  o2_mat = -o2_mat;

  std::cerr << "Computing o2 eigen decomposition...\n";
  arma::cx_vec o2_eigval_c;
  arma::vec o2_eigval(o2_dim);
  arma::cx_mat o2_eigvec_c;
  arma::mat o2_tm(o2_dim, o2_dim); // Transition matrix
  arma::eig_gen(o2_eigval_c, o2_eigvec_c, o2_mat);

  std::cerr << "o2_eigvec_c:\n" << o2_eigvec_c << std::endl;
 
  for(size_t i=0; i<o2_dim; ++i) {
    o2_eigval(i) = o2_eigval_c(i).real();
    for(size_t j=0; j<o2_dim; ++j)
      o2_tm(i,j) = o2_eigvec_c(i,j).real();
  }

  std::cerr << "o2_eigenvalues:\n"
            << o2_eigval << std::endl;
            // << "o2_eigenvectors:\n"
            // << o2_tm << std::endl;

  std::cerr << "Inverting o2 transition matrix...\n";
  arma::mat o2_inv_tm = o2_tm.i(); // Transition matrix
  std::cerr << "o2 transition matrix inverted\n";
  std::cerr << "Computing o2_tm eigen decomposition...\n";
  arma::eig_gen(o2_eigval_c, o2_eigvec_c, o2_tm);
  std::cerr << "o2_eigval_c:\n" << o2_eigval_c << std::endl;
  
  for(auto& o2_var_name : *p_o2_var_names)
    std::cout << o2_var_name << ',';
  std::cout << std::endl;

  // main loop
  std::cout << "Main loop...\n";
  std::vector<double> k_sum(o1_dim); // to store precomputed k-sum
  for(auto& t : times) {
    std::cout << "t=" << t << '\n';
    // o1
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k) {
        double l_sum=0;
        for(size_t l=0; l<o1_dim; ++l)
          l_sum += inv_tm(k,l)*c_vec(l);
        k_sum += exp(-eigval(k)*t) * tm(i,k) * l_sum;
      }
      expectations(i) = stationary_expectations(i) + k_sum;
    }
    // o2
    for(size_t i=0; i<o2_dim; ++i) {
      (*p_covariances)[i] = 0;
      for(size_t j=0; j<o2_dim; ++j) {
        double s_sum=0;
        for(size_t s=0; s<o2_dim; ++s)
          s_sum += o2_inv_tm(j,s) * (*initial_G2)(s);

        for(size_t zeta=0; zeta<o1_dim; ++zeta) { // Precomputing k-sum for all zeta (eta)
          k_sum[zeta] = 0;
          for(size_t k=0; k<o2_dim; ++k)
            k_sum[zeta] += o2_inv_tm(j,k)*(*o2_nonstationary_RHS_mat)(k,zeta);
        }

        double eta_sum=0;
        for(size_t eta=0; eta<o1_dim; ++eta)
          eta_sum += k_sum[eta]*stationary_expectations(eta);

        double mu_sum=0;
        for(size_t mu=0; mu<o1_dim; ++mu) {
          double zeta_sum=0;
          for(size_t zeta=0; zeta<o1_dim; ++zeta)
            zeta_sum += k_sum[zeta]*tm(zeta,mu);
          double nu_sum=0;
          for(size_t nu=0; nu<o1_dim; ++nu)
            nu_sum += inv_tm(mu,nu)*c_vec(nu);
          
          if(o2_eigval(j) != eigval(mu))
            mu_sum += (exp(-eigval(mu)*t)-exp(-o2_eigval(j)*t))/(o2_eigval(j)-eigval(mu)) * zeta_sum * nu_sum;
          else
            mu_sum += exp(-o2_eigval(j)*t) * zeta_sum * nu_sum;
            
          // else {
          //   std::cerr << "j=" << j << ", mu=" << mu << std::endl
          //             << "(*p_o2_var_names)[j] = " << (*p_o2_var_names)[j] << ", o1_var_names[mu] = " << o1_var_names[mu] << std::endl;            
          // }
        }
        (*p_covariances)[i] += o2_tm(i,j) * (exp(-o2_eigval(j)*t)*s_sum + (1-exp(-o2_eigval(j)*t))/o2_eigval(j)*eta_sum + mu_sum);
        // std::cerr << "t=" << t << std::endl;
        // std::cerr << "mu_sum[" << i << ',' << j << "] = " << (1-exp(-o2_eigval(j)*t))/o2_eigval(j)*eta_sum + mu_sum << std::endl; 
      }
      std::cerr << "t=" << t << ", " << (*p_o2_var_names)[i] + "=" << (*p_covariances)[i] << std::endl;
    }

    //////////// COMPUTING VARIANCES //////////
    std::vector<double> rmss(o1_dim);
    for(size_t i=0; i<o1_dim; ++i) {
      rmss[i] = ((*p_covariances)(sem_o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
      if(rmss[i]>=0) {
        std::cout << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i] << std::endl;
      }
      else
        std::cout << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i] << " NEGATIVE!\n";
    }
    // std::cout << t << ", " << expectations(2) << ", " << rmss[2] << std::endl;
    //////////////////////////////////////////
  }

  return internalise_expectations();
}


size_t Analytic_engine::o2_ind(const size_t &i, const size_t &j, const size_t &dim) const {
  if(i<=j)
    return (2*dim-i-1)*i/2+j;
  else
    return o2_ind(j, i, dim);
}

size_t Analytic_engine::sem_o2_ind(const size_t &i, const size_t &j) const {
  if(i<=j) {
    size_t n_dend = p_neuron->p_dend_segments.size(),
      n_p = 1 + n_dend + p_neuron->p_synapses.size();
    if(i==0)
      return j;
    else if(i>0 && i<n_dend+2 && j<n_dend+2)
      return o1_dim + (2+2*n_dend-i)*(i-1)/2 + j-1;
    else if(i>0 && i<n_dend+2 && j>n_dend+1)
      return (n_dend+2)*(n_dend+1)/2 + n_p*i + j;
    else
      return n_p*(n_dend+2) + (n_dend+2)*(n_dend+1)/2 + (2*n_p+n_dend-i+1)*(i-n_dend-2)/2 + j;
    }
  else
    return sem_o2_ind(j, i);
}

void Analytic_engine::initialise_o2() {
  p_o2_mat = new arma::mat(o2_dim, o2_dim);
  p_o2_RHS =  new arma::vec(o2_dim);
  p_covariances = new arma::vec(o2_dim);
  
  p_o2_var_names = new std::vector<std::string>(o2_dim);
  std::vector<std::string> &o2_var_names = *p_o2_var_names;
  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=i; j<o1_dim; ++j)
      o2_var_names[o2_ind(i,j)] = o1_var_names[i] + '-' + o1_var_names[j];
}

void Analytic_engine::set_prot_index_from(Compartment& compartment) {
  if(compartment.type() == SOMA)
    compartment.prot_ind = (p_neuron->prot_ind = p_neuron->p_dend_segments.size()+2)++;
  else
    compartment.prot_ind = p_neuron->prot_ind++;

  for (auto& p_d_comp : compartment.p_descendants)  
    set_prot_index_from(*p_d_comp);
}

void Analytic_engine::sem_initialise_o2() {
  set_prot_index_from(*p_neuron->p_soma);
  
  p_o2_mat = new arma::mat(o2_dim, o2_dim);
  p_o2_RHS =  new arma::vec(o2_dim);
  p_covariances = new arma::vec(o2_dim);
  p_o2_var_names = new std::vector<std::string>(o2_dim);

  sem_set_expectations(*sem_set_soma());

  std::vector<std::string> &o2_var_names = *p_o2_var_names;
  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=i; j<o1_dim; ++j)
      o2_var_names[sem_o2_ind(i,j)] = o1_var_names[i] + '-' + o1_var_names[j];
}

const Compartment* Analytic_engine::sem_set_soma() {
  auto& soma = *p_neuron->p_soma;

  o1_var_names[0] = soma.name + "__Gene";
  o1_var_names[1] = soma.name + "__mRNA";
  o1_var_names[soma.prot_ind] = soma.name + "__Prot";

  expectations(0) = soma.n_active_genes_expectation;
  expectations(1) = soma.n_mRNA_expectation;
  expectations(soma.prot_ind) = soma.n_prot_expectation;
  
  return &soma;
}

void Analytic_engine::sem_set_expectations(const Compartment& parent) {

  if(parent.it_p_out_junctions.empty()) {
    return;
  }

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& desc = *(*it_p_junc)->p_to;
        
    if(desc.type() == SYNAPSE) {
      o1_var_names[desc.prot_ind] = desc.name + "__Prot";

      expectations(desc.prot_ind) = desc.n_prot_expectation;
    }
    else if(desc.type() == APICAL_DENDRITE || desc.type() == BASAL_DENDRITE) {
      o1_var_names[desc.mRNA_ind] = desc.name + "__mRNA";
      o1_var_names[desc.prot_ind] = desc.name + "__Prot";

      expectations(desc.mRNA_ind) = desc.n_mRNA_expectation;
      expectations(desc.prot_ind) = desc.n_prot_expectation;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    sem_set_expectations(*((*it_p_junc)->p_to));
}

void Analytic_engine::clear_o2_matrix_and_RHS() {
  p_covariances = new arma::vec(o2_dim);
  p_o2_var_names = new std::vector<std::string>(o2_dim);
  p_o2_mat = new arma::mat(o2_dim, o2_dim);
}

void Analytic_engine::set_o2_nonstationary_RHS_soma() {

  o2_nonstationary_RHS_mat = new arma::mat(o2_dim, o1_dim);
  
  const auto& soma = *p_neuron->p_soma;

  (*o2_nonstationary_RHS_mat)(o2_ind(0,0), 0) -= soma.gene_activation_rate;
  (*o2_nonstationary_RHS_mat)(o2_ind(0,1), 0) += soma.transcription_rate;
  (*o2_nonstationary_RHS_mat)(o2_ind(1,2), 1) += soma.translation_rate;
  
  for(size_t i=0; i<o1_dim; ++i) 
    (*o2_nonstationary_RHS_mat)(o2_ind(0,i), i) += soma.gene_activation_rate*soma.number_of_gene_copies;
}

void Analytic_engine::sem_set_o2_nonstationary_RHS_soma() {

  o2_nonstationary_RHS_mat = new arma::mat(o2_dim, o1_dim);
  
  const auto& soma = *p_neuron->p_soma;

  (*o2_nonstationary_RHS_mat)(sem_o2_ind(0,0), 0) -= soma.gene_activation_rate;
  (*o2_nonstationary_RHS_mat)(sem_o2_ind(0,1), 0) += soma.transcription_rate;
  (*o2_nonstationary_RHS_mat)(sem_o2_ind(1,soma.prot_ind), 1) += soma.translation_rate;
  
  for(size_t i=0; i<o1_dim; ++i) 
    (*o2_nonstationary_RHS_mat)(sem_o2_ind(0,i), i) += soma.gene_activation_rate*soma.number_of_gene_copies;
}

void Analytic_engine::set_o2_soma() {
  const auto& soma = *p_neuron->p_soma;
  auto& o2_mat = *p_o2_mat;
  
  for(size_t i=0; i<o1_dim; ++i) {
    o2_mat(o2_ind(0,i), o2_ind(0,i)) -= (soma.gene_activation_rate + soma.gene_deactivation_rate);
    o2_mat(o2_ind(1,i), o2_ind(0,i)) += soma.transcription_rate;
    o2_mat(o2_ind(1,i), o2_ind(1,i)) -= soma.mRNA_decay_rate;
    o2_mat(o2_ind(2,i), o2_ind(1,i)) += soma.translation_rate;
    o2_mat(o2_ind(2,i), o2_ind(2,i)) -= soma.protein_decay_rate;
  }
}

void Analytic_engine::set_o2_nonstationary_RHS_mat() {
  set_o2_nonstationary_RHS_soma();

  for(auto& p_junc : p_neuron->p_junctions) { // Looping on junctions, setting descendants

    size_t& desc_start_ind = p_junc -> p_to -> o1_index;

    if(p_junc->type() != DEN_SYN)
      (*o2_nonstationary_RHS_mat)(o2_ind(desc_start_ind, desc_start_ind+1), desc_start_ind) += p_junc->p_to->translation_rate;
  }
}

void Analytic_engine::sem_set_o2_nonstationary_RHS_mat() {
  sem_set_o2_nonstationary_RHS_soma();

  for(auto& p_junc : p_neuron->p_junctions) { // Looping on junctions, setting descendants

    size_t& desc_mRNA_ind = p_junc -> p_to -> mRNA_ind;
    size_t& desc_prot_ind = p_junc -> p_to -> prot_ind;

    if(p_junc->type() != DEN_SYN)
      (*o2_nonstationary_RHS_mat)(sem_o2_ind(desc_mRNA_ind, desc_prot_ind), desc_mRNA_ind) += p_junc->p_to->translation_rate;
  }
}

void Analytic_engine::sem_set_o2_soma() {
  const auto& soma = *p_neuron->p_soma;
  auto& o2_mat = *p_o2_mat;
  auto& prot_ind = p_neuron->p_soma->prot_ind;
  
  for(size_t i=0; i<o1_dim; ++i) {
    o2_mat(sem_o2_ind(0,i), sem_o2_ind(0,i)) -= (soma.gene_activation_rate + soma.gene_deactivation_rate);
    o2_mat(sem_o2_ind(1,i), sem_o2_ind(0,i)) += soma.transcription_rate;
    o2_mat(sem_o2_ind(1,i), sem_o2_ind(1,i)) -= soma.mRNA_decay_rate;
    o2_mat(sem_o2_ind(prot_ind,i), sem_o2_ind(1,i)) += soma.translation_rate;
    o2_mat(sem_o2_ind(prot_ind,i), sem_o2_ind(prot_ind,i)) -= soma.protein_decay_rate;
  }
}

void Analytic_engine::set_o2_matrix() {
  set_o2_soma();

  auto& o2_mat = *p_o2_mat;
  
  for(auto& p_junc : p_neuron->p_junctions) { // Looping on junctions, setting descendants

    size_t& parent_start_ind = p_junc -> p_from -> o1_index;
    size_t& desc_start_ind = p_junc -> p_to -> o1_index;

    if(p_junc->type() == DEN_SYN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(parent_start_ind+1, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(parent_start_ind+1, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(desc_start_ind, i)) += p_junc->bkwd_prot_hop_rate;
      }
    else if(p_junc->type() == DEN_DEN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(o2_ind(parent_start_ind, i), o2_ind(parent_start_ind, i)) -= p_junc->fwd_mRNA_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(parent_start_ind, i)) += p_junc->fwd_mRNA_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->bkwd_mRNA_hop_rate;
        o2_mat(o2_ind(parent_start_ind, i), o2_ind(desc_start_ind, i)) += p_junc->bkwd_mRNA_hop_rate;       
        
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(parent_start_ind+1, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(parent_start_ind+1, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind+1, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(desc_start_ind+1, i)) += p_junc->bkwd_prot_hop_rate;

        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind, i)) += p_junc->p_to->translation_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind+1, i)) -= p_junc->p_to->protein_decay_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->p_to->mRNA_decay_rate;
      }
    else if(p_junc->type() == SOM_DEN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(parent_start_ind+1, i)) -= p_junc->fwd_mRNA_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(parent_start_ind+1, i)) += p_junc->fwd_mRNA_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->bkwd_mRNA_hop_rate;
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(desc_start_ind, i)) += p_junc->bkwd_mRNA_hop_rate;       
        
        o2_mat(o2_ind(parent_start_ind+2, i), o2_ind(parent_start_ind+2, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(parent_start_ind+2, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind+1, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(o2_ind(parent_start_ind+2, i), o2_ind(desc_start_ind+1, i)) += p_junc->bkwd_prot_hop_rate;

        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind, i)) += p_junc->p_to->translation_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind+1, i)) -= p_junc->p_to->protein_decay_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->p_to->mRNA_decay_rate;
      }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  auto& o2_var_names = *p_o2_var_names;
  std::cout << "-,";
  for(size_t i=0; i<o2_dim; ++i)
    std::cout << o2_var_names[i] << ',';

  for(size_t i=0; i<o2_dim; ++i) {
    std::cout << '\n';
    for(size_t j=0; j<o2_dim; ++j)
      if(j==0)
        std::cout << o2_var_names[i] << ',' << o2_mat(i,j) << ',';
      else
        std::cout << o2_mat(i,j) << ',';
  }
  std::cout << std::endl;
  // std::cout << "o2_matrix:\n" << o2_mat << std::endl;
}  


void Analytic_engine::sem_set_o2_matrix() {
  sem_set_o2_soma();

  auto& o2_mat = *p_o2_mat;
  
  for(auto& p_junc : p_neuron->p_junctions) { // Looping on junctions, setting descendants

    size_t& par_mRNA_ind = p_junc -> p_from -> mRNA_ind;
    size_t& par_prot_ind = p_junc -> p_from -> prot_ind;
    size_t& desc_mRNA_ind = p_junc -> p_to -> mRNA_ind;
    size_t& desc_prot_ind = p_junc -> p_to -> prot_ind;

    if(p_junc->type() == DEN_SYN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(par_prot_ind, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(par_prot_ind, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) += p_junc->bkwd_prot_hop_rate;
      }
    else if(p_junc->type() == DEN_DEN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(sem_o2_ind(par_mRNA_ind, i), sem_o2_ind(par_mRNA_ind, i)) -= p_junc->fwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(par_mRNA_ind, i)) += p_junc->fwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) -= p_junc->bkwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(par_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) += p_junc->bkwd_mRNA_hop_rate;       
        
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(par_prot_ind, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(par_prot_ind, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) += p_junc->bkwd_prot_hop_rate;

        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_mRNA_ind, i)) += p_junc->p_to->translation_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->p_to->protein_decay_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) -= p_junc->p_to->mRNA_decay_rate;
      }
    else if(p_junc->type() == SOM_DEN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(sem_o2_ind(par_mRNA_ind, i), sem_o2_ind(par_mRNA_ind, i)) -= p_junc->fwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(par_mRNA_ind, i)) += p_junc->fwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) -= p_junc->bkwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(par_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) += p_junc->bkwd_mRNA_hop_rate;       
        
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(par_prot_ind, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(par_prot_ind, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) += p_junc->bkwd_prot_hop_rate;

        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_mRNA_ind, i)) += p_junc->p_to->translation_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->p_to->protein_decay_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) -= p_junc->p_to->mRNA_decay_rate;
      }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }
}  


void Analytic_engine::set_o2_RHS() {
  auto& o2_RHS = *p_o2_RHS;
  const auto& soma = *p_neuron->p_soma;

  o2_RHS(o2_ind(0,0)) += soma.gene_activation_rate*expectations(0);
  o2_RHS(o2_ind(0,1)) -= soma.transcription_rate*expectations(0);
  o2_RHS(o2_ind(1,2)) -= soma.translation_rate*expectations(1);
  for(size_t i=0; i<o1_dim; ++i)
    o2_RHS(o2_ind(0,i)) -= soma.gene_activation_rate*soma.number_of_gene_copies*expectations(i);

  for(auto& ds : p_neuron->p_dend_segments)
    o2_RHS(o2_ind(ds->o1_index,ds->o1_index+1)) -= ds->translation_rate*expectations(ds->o1_index);
}

void Analytic_engine::sem_set_o2_RHS() {
  auto& o2_RHS = *p_o2_RHS;
  const auto& soma = *p_neuron->p_soma;

  o2_RHS(sem_o2_ind(0,0)) += soma.gene_activation_rate*expectations(0);
  o2_RHS(sem_o2_ind(0,1)) -= soma.transcription_rate*expectations(0);
  o2_RHS(sem_o2_ind(1,soma.prot_ind)) -= soma.translation_rate*expectations(1);
  for(size_t i=0; i<o1_dim; ++i)
    o2_RHS(sem_o2_ind(0,i)) -= soma.gene_activation_rate*soma.number_of_gene_copies*expectations(i);

  for(auto& ds : p_neuron->p_dend_segments)
    o2_RHS(sem_o2_ind(ds->mRNA_ind,ds->prot_ind)) -= ds->translation_rate*expectations(ds->mRNA_ind);
}

const Compartment* Analytic_engine::set_o2_gene_mRNA_soma() {
  Soma& soma = *(p_neuron->p_soma);
  size_t sz = 1+p_neuron->p_dend_segments.size();
  
  o2_gene_mRNA_RHS = new arma::vec(sz);
  o2_gene_mRNA_mat = new arma::mat(sz,sz);

  (*o2_gene_mRNA_RHS)(0) -= soma.gene_activation_rate*soma.number_of_gene_copies*soma.n_mRNA_expectation + soma.transcription_rate*(soma.n_active_genes_expectation + o2_gene_gene);

  (*o2_gene_mRNA_mat)(0,0) -= soma.gene_activation_rate + soma.gene_deactivation_rate + soma.mRNA_decay_rate;

  return &soma;
}

void Analytic_engine::set_o2_gene_mRNA_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty()) {
    return;
  }

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t parent_ind = p_junc -> p_from -> mRNA_ind-1;
    size_t desc_ind = p_junc -> p_to -> mRNA_ind-1;
    
    if(p_junc->type() != DEN_SYN) {

      (*o2_gene_mRNA_mat)(desc_ind, desc_ind) -= p_neuron->p_soma->gene_activation_rate + p_neuron->p_soma->gene_deactivation_rate + p_junc->p_to->mRNA_decay_rate;

      (*o2_gene_mRNA_mat)(parent_ind, parent_ind) -= p_junc->fwd_mRNA_hop_rate;
      (*o2_gene_mRNA_mat)(desc_ind, parent_ind) += p_junc->fwd_mRNA_hop_rate;
      (*o2_gene_mRNA_mat)(desc_ind, desc_ind) -= p_junc->bkwd_mRNA_hop_rate;
      (*o2_gene_mRNA_mat)(parent_ind, desc_ind) += p_junc->bkwd_mRNA_hop_rate;

      (*o2_gene_mRNA_RHS)(desc_ind) -= p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * p_junc->p_to->n_mRNA_expectation;
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o2_gene_mRNA_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::gene_mRNA_stationary_covariances() {

  set_o2_gene_mRNA_matrix(*set_o2_gene_mRNA_soma());
  
  o2_gene_mRNA = new arma::vec((*o2_gene_mRNA_mat).i()*(*o2_gene_mRNA_RHS));

  std::cerr << "o2_gene_gene = " << o2_gene_gene << std::endl;

  std::cerr << "gene_mRNA_covariances:\n" << *o2_gene_mRNA << std::endl;

  delete o2_gene_mRNA_mat; o2_gene_mRNA_mat=NULL;
  delete o2_gene_mRNA_RHS; o2_gene_mRNA_RHS=NULL;
  
  return *this;
}

const Compartment* Analytic_engine::set_o2_gene_prot_soma() {
  Soma& soma = *(p_neuron->p_soma);
  size_t sz = 1+p_neuron->p_dend_segments.size()+p_neuron->p_synapses.size();
  
  o2_gene_prot_RHS = new arma::vec(sz);
  o2_gene_prot_mat = new arma::mat(sz,sz);

  (*o2_gene_prot_RHS)(0) -= soma.gene_activation_rate*soma.number_of_gene_copies*soma.n_prot_expectation + soma.translation_rate*(*o2_gene_mRNA)(0);

  (*o2_gene_prot_mat)(0,0) -= soma.gene_activation_rate + soma.gene_deactivation_rate + soma.protein_decay_rate;

  return &soma;
}

void Analytic_engine::set_o2_gene_prot_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty()) {
    return;
  }

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t parent_ind = p_junc -> p_from -> id;
    size_t desc_ind = p_junc -> p_to -> id;
    

    (*o2_gene_prot_mat)(desc_ind, desc_ind) -= p_neuron->p_soma->gene_activation_rate + p_neuron->p_soma->gene_deactivation_rate + p_junc->p_to->protein_decay_rate;

    (*o2_gene_prot_mat)(parent_ind, parent_ind) -= p_junc->fwd_prot_hop_rate;
    (*o2_gene_prot_mat)(desc_ind, parent_ind) += p_junc->fwd_prot_hop_rate;
    (*o2_gene_prot_mat)(desc_ind, desc_ind) -= p_junc->bkwd_prot_hop_rate;
    (*o2_gene_prot_mat)(parent_ind, desc_ind) += p_junc->bkwd_prot_hop_rate;
    
    if(p_junc->type() != DEN_SYN)
      (*o2_gene_prot_RHS)(desc_ind) -= p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * p_junc->p_to->n_prot_expectation + p_junc->p_to->translation_rate*(*o2_gene_mRNA)(p_junc->p_to->mRNA_ind-1);
    else
      (*o2_gene_prot_RHS)(desc_ind) -= p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * p_junc->p_to->n_prot_expectation;
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o2_gene_prot_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::gene_protein_stationary_covariances() {

  set_o2_gene_prot_matrix(*set_o2_gene_prot_soma());
  
  o2_gene_prot = new arma::vec((*o2_gene_prot_mat).i()*(*o2_gene_prot_RHS));

  std::cerr << "gene_prot_covariances:\n" << *o2_gene_prot << std::endl;

  delete o2_gene_prot_mat; o2_gene_prot_mat=NULL;
  delete o2_gene_prot_RHS; o2_gene_prot_RHS=NULL;

  return *this;
}

const Compartment* Analytic_engine::set_o2_mRNA_mRNA_soma() {
  Soma& soma = *(p_neuron->p_soma);
  size_t sz = 1+p_neuron->p_dend_segments.size();
  
  o2_mRNA_mRNA_RHS = new arma::vec(sz*(sz+1)/2);
  o2_mRNA_mRNA_mat = new arma::mat(sz*(sz+1)/2, sz*(sz+1)/2);

  for(size_t i=0; i<sz; ++i) {
    (*o2_mRNA_mRNA_RHS)(i) -= soma.transcription_rate*(*o2_gene_mRNA)(i);
    (*o2_mRNA_mRNA_mat)(i,i) -= soma.mRNA_decay_rate;
  }

  return &soma;
}

void Analytic_engine::set_o2_mRNA_mRNA_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty()) {
    return;
  }

  size_t sz = 1 + p_neuron->p_dend_segments.size();

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;

    if(p_junc->type() != DEN_SYN) {
      
      size_t parent_ind = p_junc -> p_from -> mRNA_ind-1;
      size_t desc_ind = p_junc -> p_to -> mRNA_ind-1;
    
      for(size_t i=0; i<sz; ++i) {
        (*o2_mRNA_mRNA_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->p_to->mRNA_decay_rate;

        (*o2_mRNA_mRNA_mat)(o2_ind(parent_ind, i, sz), o2_ind(parent_ind, i, sz)) -= p_junc->fwd_mRNA_hop_rate;
        (*o2_mRNA_mRNA_mat)(o2_ind(desc_ind, i, sz), o2_ind(parent_ind, i, sz)) += p_junc->fwd_mRNA_hop_rate;
        (*o2_mRNA_mRNA_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->bkwd_mRNA_hop_rate;
        (*o2_mRNA_mRNA_mat)(o2_ind(parent_ind, i, sz), o2_ind(desc_ind, i, sz)) += p_junc->bkwd_mRNA_hop_rate;
      }
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o2_mRNA_mRNA_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::mRNA_mRNA_stationary_covariances() {
  std::cout << "Computing mRNA_mRNA_stationary_covariances...\n";
  set_o2_mRNA_mRNA_matrix(*set_o2_mRNA_mRNA_soma());
  
  o2_mRNA_mRNA = new arma::vec((*o2_mRNA_mRNA_mat).i()*(*o2_mRNA_mRNA_RHS));

  std::cerr << "mRNA_mRNA_covariances:\n" << *o2_mRNA_mRNA << std::endl;

  delete o2_mRNA_mRNA_mat; o2_mRNA_mRNA_mat=NULL;
  delete o2_mRNA_mRNA_RHS; o2_mRNA_mRNA_RHS=NULL;

  return *this;
}

const Compartment* Analytic_engine::set_o2_mRNA_prot_soma() {
  Soma& soma = *(p_neuron->p_soma);
  size_t sz_mRNA = 1+p_neuron->p_dend_segments.size();
  size_t sz_prot = sz_mRNA + p_neuron->p_synapses.size();
  
  o2_mRNA_prot_RHS = new arma::vec(sz_mRNA*sz_prot);
  o2_mRNA_prot_mat = new arma::mat(sz_mRNA*sz_prot, sz_mRNA*sz_prot);

  (*o2_mRNA_prot_RHS)(0) -= soma.translation_rate*soma.n_mRNA_expectation; // 0==o2_ind_asym(0,0,sz_prot)

  for(size_t i=0; i<sz_mRNA; ++i) {
    (*o2_mRNA_prot_RHS)(i*sz_prot) -= soma.translation_rate*(*o2_mRNA_mRNA)(i); //i==o2_ind(0,i,sz_mRNA)
    (*o2_mRNA_prot_mat)(i*sz_prot, i*sz_prot) -= soma.protein_decay_rate; // i*sz_prot==o2_ind_asym(i,0,sz_mRNA)
  }
  for(size_t i=0; i<sz_prot; ++i) {
    (*o2_mRNA_prot_RHS)(i) -= soma.transcription_rate*(*o2_gene_prot)(i);
    (*o2_mRNA_prot_mat)(i, i) -= soma.mRNA_decay_rate; // i==o2_ind_asym(0,i,sz_prot)
  }
  return &soma;
}

void Analytic_engine::set_o2_mRNA_prot_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty()) {
    return;
  }

  size_t sz_mRNA = 1 + p_neuron->p_dend_segments.size(),
    sz_prot = sz_mRNA + p_neuron->p_synapses.size();
  
  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;

    size_t parent_mRNA_ind = p_junc -> p_from -> mRNA_ind-1,
      desc_mRNA_ind = p_junc -> p_to -> mRNA_ind-1,
      parent_prot_ind = p_junc -> p_from -> id,
      desc_prot_ind = p_junc -> p_to -> id;

    if(p_junc->type() != DEN_SYN) {

      (*o2_mRNA_prot_RHS)(o2_ind_asym(desc_mRNA_ind, desc_prot_ind, sz_prot)) -= (p_junc->p_to->translation_rate)*(p_junc->p_to->n_mRNA_expectation);
    
      for(size_t i=0; i<sz_mRNA; ++i) {
        (*o2_mRNA_prot_RHS)(o2_ind_asym(i, desc_prot_ind, sz_prot)) -= (p_junc->p_to->translation_rate)*(*o2_mRNA_mRNA)(o2_ind(desc_mRNA_ind, i, sz_mRNA));
        
        (*o2_mRNA_prot_mat)(o2_ind_asym(i,desc_prot_ind,sz_prot), o2_ind_asym(i,desc_prot_ind,sz_prot)) -= p_junc->p_to->protein_decay_rate;

        (*o2_mRNA_prot_mat)(o2_ind_asym(i, parent_prot_ind, sz_prot), o2_ind_asym(i, parent_prot_ind, sz_prot)) -= p_junc->fwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, desc_prot_ind, sz_prot), o2_ind_asym(i, parent_prot_ind, sz_prot)) += p_junc->fwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, desc_prot_ind, sz_prot), o2_ind_asym(i, desc_prot_ind, sz_prot)) -= p_junc->bkwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, parent_prot_ind, sz_prot), o2_ind_asym(i, desc_prot_ind, sz_prot)) += p_junc->bkwd_prot_hop_rate;
      }
      for(size_t i=0; i<sz_prot; ++i) {
        (*o2_mRNA_prot_mat)(o2_ind_asym(desc_mRNA_ind,i,sz_prot), o2_ind_asym(desc_mRNA_ind,i,sz_prot)) -= p_junc->p_to->mRNA_decay_rate;

        (*o2_mRNA_prot_mat)(o2_ind_asym(parent_mRNA_ind, i, sz_prot), o2_ind_asym(parent_mRNA_ind, i, sz_prot)) -= p_junc->fwd_mRNA_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(desc_mRNA_ind, i, sz_prot), o2_ind_asym(parent_mRNA_ind, i, sz_prot)) += p_junc->fwd_mRNA_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(desc_mRNA_ind, i, sz_prot), o2_ind_asym(desc_mRNA_ind, i, sz_prot)) -= p_junc->bkwd_mRNA_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(parent_mRNA_ind, i, sz_prot), o2_ind_asym(desc_mRNA_ind, i, sz_prot)) += p_junc->bkwd_mRNA_hop_rate;
      }
    }
    else {
      for(size_t i=0; i<sz_mRNA; ++i) {
        (*o2_mRNA_prot_mat)(o2_ind_asym(i,desc_prot_ind,sz_prot), o2_ind_asym(i,desc_prot_ind,sz_prot)) -= p_junc->p_to->protein_decay_rate;

        (*o2_mRNA_prot_mat)(o2_ind_asym(i, parent_prot_ind, sz_prot), o2_ind_asym(i, parent_prot_ind, sz_prot)) -= p_junc->fwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, desc_prot_ind, sz_prot), o2_ind_asym(i, parent_prot_ind, sz_prot)) += p_junc->fwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, desc_prot_ind, sz_prot), o2_ind_asym(i, desc_prot_ind, sz_prot)) -= p_junc->bkwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, parent_prot_ind, sz_prot), o2_ind_asym(i, desc_prot_ind, sz_prot)) += p_junc->bkwd_prot_hop_rate;
      }
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o2_mRNA_prot_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::mRNA_protein_stationary_covariances() {
  std::cout << "Computing mRNA_protein_stationary_covariances...\n";
  set_o2_mRNA_prot_matrix(*set_o2_mRNA_prot_soma());
  
  o2_mRNA_prot = new arma::vec((*o2_mRNA_prot_mat).i()*(*o2_mRNA_prot_RHS));

  std::cerr << "mRNA_prot_covariances:\n" << *o2_mRNA_prot << std::endl;

  delete o2_mRNA_prot_mat; o2_mRNA_prot_mat=NULL;
  delete o2_mRNA_prot_RHS; o2_mRNA_prot_RHS=NULL;

  return *this;
}

const Compartment* Analytic_engine::set_o2_prot_prot_soma() {
  Soma& soma = *(p_neuron->p_soma);
  size_t sz = 1 + p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size();
  
  o2_prot_prot_RHS = new arma::vec(sz*(sz+1)/2);
  o2_prot_prot_mat = new arma::mat(sz*(sz+1)/2, sz*(sz+1)/2);

  for(size_t i=0; i<sz; ++i) {
    (*o2_prot_prot_RHS)(i) -= soma.translation_rate*(*o2_mRNA_prot)(i);
    (*o2_prot_prot_mat)(i,i) -= soma.protein_decay_rate;
  }

  return &soma;
}

void Analytic_engine::set_o2_prot_prot_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty()) {
    return;
  }

  size_t sz = 1 + p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size();

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;

    size_t parent_ind = p_junc -> p_from -> id;
    size_t desc_ind = p_junc -> p_to -> id;
    
    if(p_junc->type() != DEN_SYN) {
      size_t desc_mRNA_ind = p_junc -> p_to -> mRNA_ind-1;
    
      for(size_t i=0; i<sz; ++i) {
        (*o2_prot_prot_RHS)(o2_ind(desc_ind, i, sz)) -= p_junc->p_to->translation_rate*(*o2_mRNA_prot)(o2_ind_asym(desc_mRNA_ind, i, sz));
        
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->p_to->protein_decay_rate;

        (*o2_prot_prot_mat)(o2_ind(parent_ind, i, sz), o2_ind(parent_ind, i, sz)) -= p_junc->fwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(parent_ind, i, sz)) += p_junc->fwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->bkwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(parent_ind, i, sz), o2_ind(desc_ind, i, sz)) += p_junc->bkwd_prot_hop_rate;
      }
    }
    else {
      for(size_t i=0; i<sz; ++i) {        
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->p_to->protein_decay_rate;

        (*o2_prot_prot_mat)(o2_ind(parent_ind, i, sz), o2_ind(parent_ind, i, sz)) -= p_junc->fwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(parent_ind, i, sz)) += p_junc->fwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->bkwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(parent_ind, i, sz), o2_ind(desc_ind, i, sz)) += p_junc->bkwd_prot_hop_rate;
      }      
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o2_prot_prot_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::protein_protein_stationary_covariances() {
  std::cout << "Computing protein_protein_stationary_covariances...\n";
  set_o2_prot_prot_matrix(*set_o2_prot_prot_soma());
  
  o2_prot_prot = new arma::vec((*o2_prot_prot_mat).i()*(*o2_prot_prot_RHS));

  std::cerr << "protein_protein_covariances:\n" << *o2_prot_prot << std::endl;

  delete o2_prot_prot_mat; o2_prot_prot_mat=NULL;
  delete o2_prot_prot_RHS; o2_prot_prot_RHS=NULL;
  
  return *this;
}

Analytic_engine& Analytic_engine::stationary_covariances() {
  initialise_o2();
  set_o2_matrix();
  set_o2_RHS();
  std::cerr << "(*p_o2_RHS)(j) = " << (*p_o2_RHS) << std::endl;
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  auto& o2_RHS = *p_o2_RHS;
  std::cerr << "Computing covariances...\n";
  covariances = o2_mat.i()*o2_RHS;

  std::vector<double> rmss(o1_dim);
  for(size_t i=0; i<o1_dim; ++i) {  
    rmss[i] = sqrt(covariances(o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
    std::cerr << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i] << std::endl;
  }
  
  return *this;
}

Analytic_engine& Analytic_engine::sem_stationary_covariances() {
  sem_initialise_o2();
  sem_set_o2_matrix();
  sem_set_o2_RHS();
  
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  auto& o2_RHS = *p_o2_RHS;
  std::cerr << "Computing covariances...\n";
  covariances = o2_mat.i()*o2_RHS;

  std::vector<double> rmss(o1_dim);
  for(size_t i=0; i<o1_dim; ++i) {  
    rmss[i] = sqrt(covariances(sem_o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
    std::cerr << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i]<< ", " << rmss[i]/expectations(i) << std::endl;
  }

  std::cerr << "--,";
  for(unsigned int i=0; i<o1_dim; ++i)
    std::cerr << o1_var_names[i] << ',';
  for(unsigned int i=0; i<o1_dim; ++i) {
    std::cerr << std::endl << o1_var_names[i] << ',';
    for(unsigned int j=0; j<o1_dim; ++j)
      std::cerr << covariances(sem_o2_ind(i,j)) << ',';
  }
  std::cerr << std::endl;
  
  return *this;
}

Analytic_engine& Analytic_engine::sem_stationary_pearson_correlations() {

  sem_stationary_covariances();
  auto& covariances = *p_covariances;

  std::vector<double> rmss(o1_dim);
  for(size_t i=0; i<o1_dim; ++i) {  
    rmss[i] = sqrt(covariances(sem_o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
    std::cerr << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i]<< ", " << rmss[i]/expectations(i) << std::endl;
  }

  std::cerr << "\n----------------------------------------------------\n";

  std::cerr << "--,";
  for(unsigned int i=0; i<o1_dim; ++i)
    std::cerr << o1_var_names[i] << ',';
  for(unsigned int i=0; i<o1_dim; ++i) {
    std::cerr << std::endl << o1_var_names[i] << ',';
    for(unsigned int j=0; j<o1_dim; ++j)
      if(i != j)
        std::cerr << (covariances(sem_o2_ind(i,j)) - expectations(i)*expectations(j))/(rmss[i]*rmss[j]) << ',';
      else
        std::cerr << "1,";
  }
  std::cerr << std::endl;
  
  return *this;
}

Analytic_engine& Analytic_engine::clear_all() {
  if(p_o2_mat) delete p_o2_mat;
  if(p_o2_RHS) delete p_o2_RHS;
  if(p_covariances) delete p_covariances;
  if(p_o2_var_names) delete p_o2_var_names;
  if(o2_gene_mRNA) delete o2_gene_mRNA;
  if(o2_gene_mRNA_RHS) delete o2_gene_mRNA_RHS;
  if(o2_gene_mRNA_mat) delete o2_gene_mRNA_mat;
  if(o2_gene_prot) delete o2_gene_prot;
  if(o2_gene_prot_RHS) delete o2_gene_prot_RHS;
  if(o2_gene_prot_mat) delete o2_gene_prot_mat;
  if(o2_mRNA_mRNA) delete o2_mRNA_mRNA;
  if(o2_mRNA_mRNA_RHS) delete o2_mRNA_mRNA_RHS;
  if(o2_mRNA_mRNA_mat) delete o2_mRNA_mRNA_mat;
  if(o2_mRNA_prot) delete o2_mRNA_prot;
  if(o2_mRNA_prot_RHS) delete o2_mRNA_prot_RHS;
  if(o2_mRNA_prot_mat) delete o2_mRNA_prot_mat;
  if(o2_nonstationary_RHS_mat) delete o2_nonstationary_RHS_mat;

  return *this;
}

Analytic_engine::~Analytic_engine() {
  clear_all();
}
