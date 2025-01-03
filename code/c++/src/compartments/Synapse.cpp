#include "../../include/compartments/Synapse.hpp"

Synapse::Synapse(Compartment &parent, const std::string& name) : Compartment(2.5, name) {

  protein_decay_rate = 0;

  parent.p_descendants.push_back(this);
  
  if(parent.p_neuron) {
    p_neuron = parent.p_neuron;
    p_neuron -> p_synapses.push_back(this);
  }
}

Synapse::Synapse(Compartment &parent, const std::string& name, const double &protein_binding_rate, const double &protein_unbinding_rate, const double &protein_decay_rate, const unsigned int &protein_number) : Compartment(2.5, name), protein_binding_rate(protein_binding_rate), protein_unbinding_rate(protein_unbinding_rate) {
  this->protein_decay_rate = protein_decay_rate;
  
  parent.p_descendants.push_back(this);
  
  if(parent.p_neuron) {
    p_neuron = parent.p_neuron;
    p_neuron -> p_synapses.push_back(this);
  }
}

std::ostream& operator<<(std::ostream& os, const Synapse& syn) {
  os << syn.get_name() << ", params: " << "br=" << syn.protein_binding_rate << "; ubr=" << syn.protein_unbinding_rate;
  return os;
}
