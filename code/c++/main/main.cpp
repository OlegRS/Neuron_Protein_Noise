#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Synapse.hpp"
#include "../engines/analytic/Analytic_engine.hpp"
#include "../engines/stochastic/Gillespie_engine.hpp"

using namespace std;

#define N_FORKS 1 // Note that it is (2^N_FORKS - 1)*3 compartments (if 2 synapses on each dend seg)!
void fork_dendrite(Dendritic_segment* ds, size_t depth=0) {
  if (depth < N_FORKS) {
    auto ds1 = new Dendritic_segment(*ds, ds->get_name() + "-1");
    new Synapse(*ds1, "s_" + ds1->get_name() + "_1");
    new Synapse(*ds1, "s_" + ds1->get_name() + "_2");
    fork_dendrite(ds1, depth+1);

    auto ds2 = new Dendritic_segment(*ds, ds->get_name() + "-2");
    new Synapse(*ds2, "s_" + ds2->get_name() + "_1");
    new Synapse(*ds2, "s_" + ds2->get_name() + "_2", .6, 6 + .00000001);
    fork_dendrite(ds2, depth+1);
  }
}

int main() {

  Soma soma("soma" /*,Parameters of the soma*/);
  
  ///// Linear neuron
  // Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1");
  // new Synapse(*p_ds, "s_1_1");
  // new Synapse(*p_ds, "s_1_2");
  
  // for(unsigned int i=0; i<100; ++i) {
  //   p_ds = new Dendritic_segment(*p_ds, "d_" + to_string(i+2));
  //   new Synapse(*p_ds, "s_" + to_string(i+2) + "_1");
  //   new Synapse(*p_ds, "s_" + to_string(i+2) + "_2");
  // }


  ///// Branching neuron
  Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1");
  new Synapse(*p_ds, "s_1_1");
  new Synapse(*p_ds, "s_1_2");

  fork_dendrite(p_ds);


  Neuron neuron(soma, "Test_neuron");
  
  // cout << neuron << endl;


  // cout << "----------------- ANALYTIC ENGINE -----------------\n";
  Analytic_engine ae(neuron);
  // ae.stationary_expectations().stationary_covariances();
  // ae.stationary_expectations().sem_stationary_covariances();
  // ae.stationary_expectations().sem_stationary_pearson_correlations();

  // ae.mRNA_stationary_expectations().protein_stationary_expectations();
  // ae.mRNA_o1_eigen_decomposition();
  // ae.protein_o1_eigen_decomposition()

  // ae.gene_mRNA_stationary_covariances().gene_protein_stationary_covariances().mRNA_mRNA_stationary_covariances().mRNA_protein_stationary_covariances().protein_protein_stationary_covariances();

  std::list<double> times;
  for(size_t i=0; i<10000; ++i)
    times.push_back(i);

  // Analytic_engine(neuron).sem_nonstationary_expectations(times);

  ae.sem_stationary_expectations();
  // ae.stationary_covariances();

  ae.sem_nonstationary_covariances(times);
  // Gillespie_engine(neuron).run_Gillespie(10000);

  // ae.sem_nonstationary_expectations(times);
  
  return 0;
}