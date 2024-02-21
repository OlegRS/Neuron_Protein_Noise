#include "../compartments/Synapse.hpp"

Synapse::Synapse(Compartment &parent, const std::string& name) : Compartment(name) {

  protein_decay_rate *= 10;

  parent.p_descendants.push_back(this);
  
  if(parent.p_neuron) {
    p_neuron = parent.p_neuron;
    p_neuron -> p_synapses.push_back(this);
  }
}

Synapse::Synapse(Compartment &parent, const std::string& name, const double &protein_binding_rate, const double &protein_unbinding_rate, const unsigned int &protein_number) : Compartment(name), protein_binding_rate(protein_binding_rate), protein_unbinding_rate(protein_unbinding_rate) {

  protein_decay_rate *= 10;
  
  parent.p_descendants.push_back(this);
  
  if(parent.p_neuron) {
    p_neuron = parent.p_neuron;
    p_neuron -> p_synapses.push_back(this);
  }
}


std::ostream& operator<<(std::ostream& os, const Synapse& syn) {
  os << syn.get_name() << ", params: ";
  return os;
}
