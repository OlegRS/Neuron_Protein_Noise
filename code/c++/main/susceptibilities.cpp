#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Synapse.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

int main() {

  Soma soma("soma" /*,Parameters of the soma*/);

  ///// Branching neuron
  std::list<Synapse*> p_synapses;
  Dendritic_segment ds(soma, "d_1");
  Synapse syn_1_1(ds, "s_1_1", .6, 6, 1.2e-5*3600 * 10);
  p_synapses.push_back(&syn_1_1);
  Synapse syn_1_2(ds, "s_1_2", .6, 6, 1.2e-5*3600 * 10);
  p_synapses.push_back(&syn_1_2);
  Dendritic_segment ds_1(ds, "d_1-1");
  Synapse syn_11_1(ds_1, "s_1_1-1", .6, 6, 1.2e-5*3600 * 10);
  p_synapses.push_back(&syn_11_1);
  Synapse syn_11_2(ds_1, "s_1_1-2", .6, 6, 1.2e-5*3600 * 10);
  p_synapses.push_back(&syn_11_2);
  Dendritic_segment ds_2(ds, "d_1-2");
  Synapse syn_12_1(ds_2, "s_1_2-1", .6, 6, 1.2e-5*3600 * 10);
  p_synapses.push_back(&syn_12_1);
  Synapse syn_12_2(ds_2, "s_1_2-2", .6, 6, 1.2e-5*3600 * 10);
  p_synapses.push_back(&syn_12_2);

  Neuron neuron(soma, "Test_neuron");
  
  std::cout << neuron << std::endl;
  
  Analytic_engine ae(neuron);

#define dim 15
  
#define pdr_start 0 // protein decay rate
#define pdr_fin   50
#define d_prot_dec_rate     .0001
#define prot_bind_rate      .6
#define d_prot_bind_rate    .001
#define SYN                 syn_11_2

  std::ofstream ofs_exp_diff("expectation_differences"),
    ofs_expectations1("expectations1"),
    ofs_expectations2("expectations2");

  arma::vec exp_diff(dim), expectations(dim);
  
  std::cerr << "------------------- Merged loop -----------------------\n";
  for(double prot_dec_rate=pdr_start; prot_dec_rate<pdr_fin; prot_dec_rate+=d_prot_dec_rate) {
    for(auto syn : p_synapses)
      syn->set_protein_decay_rate(prot_dec_rate);
    //    SYN.set_protein_decay_rate(prot_dec_rate);
    
    SYN.set_protein_binding_rate(prot_bind_rate + d_prot_bind_rate);
    neuron.refresh();
    ae.stationary_expectations();
    exp_diff  = ae.get_expectations();

    ofs_expectations1 << prot_dec_rate;
    for(size_t i=0; i<dim; ++i)
      ofs_expectations1  << ',' << exp_diff(i);
    ofs_expectations1 << std::endl;
    
 
    SYN.set_protein_binding_rate(prot_bind_rate - d_prot_bind_rate);
    neuron.refresh();
    ae.stationary_expectations();

    exp_diff -= expectations = ae.get_expectations();
    ofs_expectations2 << prot_dec_rate;
    for(size_t i=0; i<dim; ++i)
      ofs_expectations2  << ',' << expectations(i);
    ofs_expectations2 << std::endl;

    for(size_t i=0; i<dim; ++i)
      exp_diff[i] /= (2*d_prot_bind_rate) * expectations[i];
    ofs_exp_diff << prot_dec_rate;
    for(size_t i=0; i<dim; ++i)
      ofs_exp_diff  << ',' << exp_diff(i);
    ofs_exp_diff << std::endl;
  }

  auto var_names = *(ae.o1_variable_names());
  for(size_t i=0; i<dim; ++i)
    std::cerr << i+1 << ") " << var_names[i] << '\n';
  
  return 0;
}
