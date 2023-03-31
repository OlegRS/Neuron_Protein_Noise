#ifndef __NEURON_HPP__
#define __NEURON_HPP__

#include <iostream>
#include <list>

#include "compartments/Compartment.hpp"
#include "junctions/Junction.hpp"

class Soma;
class Dendritic_segment;
class Synapse;

class Analytic_engine;
class Gillespie_engine;
class Downsampling_engine;

class Neuron {
  friend class Compartment;
  friend class Soma;
  friend class Dendritic_segment;
  friend class Synapse;
  friend class Analytic_engine;
  friend class Gillespie_engine;
  friend class Downsampling_engine;

protected:

  class Som_den_junction;
  class Den_syn_junction;
  class Den_den_junction;

  // Auxiliary variables for Analytic_engine
  size_t o1_index = 0, mRNA_ind=1, prot_ind, comp_id=0;
  // For Gillespie_engine
  double total_rate = 0; // Combined rate of events
  
  std::string name;

  Soma* p_soma = NULL;
  std::list<Compartment*> p_dend_segments;
  std::list<Compartment*> p_synapses;
  
  std::list<Junction*> p_junctions;
  unsigned int n_SDJ, n_DSJ, n_DDJ; // Numbers of different junctions

  void associate(Compartment& comaprtment); // Associates compartment with the neuron
  void dissociate(Compartment& comaprtment); // Dissociates compartment from the neuron

public:
 
  Neuron() {};
  Neuron(Soma &soma, const std::string& name="no_name");
  Neuron(const Neuron&);
  Neuron(const std::string& file_name);

  void save(const std::string& file_name) const;
  Neuron& load(const std::string& file_name);

  Compartment& add_compartment(Compartment&);

  Neuron& connect(Compartment& comp1, Compartment& comp2);

  Neuron& link_compartments(Compartment &comp1, Compartment &comp2);

  friend std::ostream& operator<<(std::ostream&, const Neuron&);
  friend std::ostream& operator<<(std::ostream&, const Junction&);

  ~Neuron();
};

#endif