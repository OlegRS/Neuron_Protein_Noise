#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Synapse.hpp"
#include "../engines/stochastic/Gillespie_engine.hpp"

#define N_AVRG  10000
#define t1 5000
#define t2 10000

#define N_FORKS 1 // Note that it is (2^N_FORKS - 1)*3 compartments (if 2 synapses on each dend seg)!
void fork_dendrite(Dendritic_segment* ds, size_t depth=0) {
  if (depth < N_FORKS) {
    auto ds1 = new Dendritic_segment(*ds, ds->get_name() + "-1");
    new Synapse(*ds1, "s_" + ds1->get_name() + "_1");
    new Synapse(*ds1, "s_" + ds1->get_name() + "_2", .6, 6);
    fork_dendrite(ds1, depth+1);

    auto ds2 = new Dendritic_segment(*ds, ds->get_name() + "-2");
    new Synapse(*ds2, "s_" + ds2->get_name() + "_1");
    new Synapse(*ds2, "s_" + ds2->get_name() + "_2", .6, 6);
    fork_dendrite(ds2, depth+1);
  }
}

int main() {

  PRNG rnd(1);

  double dt = .1;  
  
  std::string file_name = "../../data/gillespie/heterosynaptic_plasticity/HM_";

  std::list<double> times1;
  for(double t=0; t<t1; t+=dt)
    times1.push_back(t);

  std::list<double> times2;
  for(double t=t1; t<t2; t+=dt)
    times2.push_back(t);

  
  for(size_t i=0; i<N_AVRG; ++i) {
    Soma soma("soma" /*,Parameters of the soma*/);
    // ///// Branching neuron
    // Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1");
    // Synapse *p_syn_1_1 = new Synapse(*p_ds, "s_1_1", .6, 6);
    // Synapse *p_syn_1_2 = new Synapse(*p_ds, "s_1_2", .6, 6);
    // // fork_dendrite(p_ds);

    Dendritic_segment ds(soma, "d_1");
    Synapse syn_1_1(ds, "s_1_1", .6, 6);
    Synapse syn_1_2(ds, "s_1_1", .6, 6);

    Neuron neuron(soma, "Test_neuron");  

    std::ofstream ofs_gillespie(file_name + std::to_string(i));

    std::cerr << "Writing Gillespie results to: " << file_name + std::to_string(i) << '\n';
  
    std::cerr << "------------------- Loop_1 -----------------------\n";
    std::cout << neuron << std::endl;
 
    Gillespie_engine(neuron, rnd).run_Gillespie(times1, ofs_gillespie);
  
    std::cerr << "------------------- Loop_2 -----------------------\n";

    // syn_1_2.set_protein_binding_rate(.6);
    // ds.set_translation_rate(0.021*3600*10);
    
    neuron.refresh();
    std::cout << neuron << std::endl;
  
    Gillespie_engine(neuron, rnd).run_Gillespie(times2, ofs_gillespie);
  
    ofs_gillespie.close();
  }
  
  return 0;
}
