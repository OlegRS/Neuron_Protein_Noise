#include "Neuron.hpp"
#include "compartments/Soma.hpp"
#include "compartments/Dendritic_segment.hpp"
#include "compartments/Synapse.hpp"
#include "junctions/Som_den_junction.hpp"
#include "junctions/Den_den_junction.hpp"
#include "junctions/Den_syn_junction.hpp"

void Neuron::associate(Compartment& compartment) {

  compartment.p_neuron = this;

  Compartment::Type p_type = compartment.type();
  compartment.id = comp_id++;
  
  if(p_type == SYNAPSE) {
    p_synapses.push_back(&compartment);
    compartment.iterator = p_synapses.end();
    compartment.o1_index = o1_index;
    o1_index += 1;
  }
  else if(p_type == APICAL_DENDRITE || p_type == BASAL_DENDRITE) {
    p_dend_segments.push_back(&compartment);
    compartment.iterator = p_dend_segments.end();
    compartment.o1_index = o1_index;
    o1_index += 2;
    compartment.mRNA_ind = mRNA_ind++;
  }
  else if(p_type == SOMA) {
    compartment.o1_index = o1_index;
    o1_index += 3;
    compartment.mRNA_ind = mRNA_ind++;
  }

  // Writing junctions
  for (auto& p_d_comp : compartment.p_descendants){
    auto d_type = p_d_comp->type();
    bool p_dend = p_type == APICAL_DENDRITE || p_type == BASAL_DENDRITE;
    bool d_dend = d_type == APICAL_DENDRITE || d_type == BASAL_DENDRITE;
    if(p_type == SOMA && d_dend) {
      p_junctions.push_back(new Som_den_junction(&compartment, p_d_comp));
      compartment.it_p_out_junctions.push_back(--p_junctions.end());
      p_d_comp->it_p_in_junctions.push_back(--p_junctions.end());
      ++n_SDJ;
    }
    else if(p_dend && d_type==SYNAPSE) {
      p_junctions.push_back(new Den_syn_junction(&compartment, p_d_comp));
      compartment.it_p_out_junctions.push_back(--p_junctions.end());
      p_d_comp->it_p_in_junctions.push_back(--p_junctions.end());
      ++n_DSJ;
    }
    else if(p_dend && d_dend) {
      p_junctions.push_back(new Den_den_junction(&compartment, p_d_comp));
      compartment.it_p_out_junctions.push_back(--p_junctions.end());
      p_d_comp->it_p_in_junctions.push_back(--p_junctions.end());
      ++n_DDJ;
    }
    else {
      std::cerr << "---------------------------------------------------\n"
                << "ERROR: " << d_type << " descending from " << p_type << std::endl
                << "---------------------------------------------------\n";
      exit(1);
    }
    
    associate(*p_d_comp);
  }
}

// void Neuron::dissociate(Compartment& compartment) {
//   compartment.p_neuron = NULL;

//   Compartment::Type type = compartment.type();

//   if(type.id == SYNAPSE) 
//     p_synapses.erase(compartment.iterator);
//   else if(type.id == APICAL_DENDRITE || type.id == BASAL_DENDRITE)
//     p_dend_segments.erase(compartment.iterator);
  
//   for (auto& comp : compartment.p_descendants)
//     dissociate(*comp);
// }

Neuron::Neuron(Soma &soma, const std::string &name) : name(name) {
  p_soma = &soma;
  associate(static_cast<Compartment&>(soma));
}

void Neuron::save(const std::string& file_name) const {
}

std::ostream& operator<<(std::ostream &os , const Neuron &neur) {

  os << "******* COMPARTMENTS *******:\n"
     << "* Soma:\n** " << *(Compartment*)neur.p_soma << std::endl
     << "\n* Dendritic segments:\n";

  for(auto& p_ds : neur.p_dend_segments)
    os << "** " << *p_ds << "\n";
  os << '\n';
  os << "* Synapses:\n";
  for(auto& p_syn : neur.p_synapses)
    os << "** " << *p_syn << "\n";
  os << '\n';
  os << "******* JUNCTIONS *******:\n";
  for(auto& junct : neur.p_junctions)
    os << "** " << *junct << "\n";
  
  return os;
}


Neuron::~Neuron() {
  for(auto& junct : p_junctions)
    delete junct;
}
