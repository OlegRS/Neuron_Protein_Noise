#include "../../include/compartments/Dendritic_segment.hpp"
#include "../../include/compartments/Soma.hpp"

Dendritic_segment::Dendritic_segment(Compartment& parent, const std::string& name, const double& length, const double& radius, const double& theta, const double& phi) : Compartment(parent, name, length, radius, theta, phi)  {
  
  parent.p_descendants.push_back(this);

  if(parent.type() == BASAL_DENDRITE || parent.type() == APICAL_DENDRITE)
    ++static_cast<Dendritic_segment&>(parent).n_descending_DS;
  else if(parent.type() == SOMA)
    ++static_cast<Soma&>(parent).n_descending_DS;
  
  if(parent.p_neuron) {
    p_neuron = parent.p_neuron;
    p_neuron -> p_dend_segments.push_back(this);
  }
}

std::ostream& operator<<(std::ostream& os, const Dendritic_segment& ds) {
  os << ds.get_name() <<", n_descending_DS: "<< ds.n_descending_DS;
  os << "\nIn_junctions:\n";
  for(auto& it_p_in_junc : ds.it_p_in_junctions)
    os << **it_p_in_junc << std::endl;
  os << "\nOut_junctions:\n";
  for(auto& it_p_out_junc : ds.it_p_out_junctions)
    os << **it_p_out_junc << std::endl;
  return os;
}
