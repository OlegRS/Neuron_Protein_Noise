#include "Compartment.hpp"

Compartment::Type::operator std::string() const {
  if(id==SOMA) return "Soma";
  else if(id==BASAL_DENDRITE) return "Basal dendritic segment";
  else if(id==APICAL_DENDRITE) return "Apical dendritic segment";
  else if(id==SYNAPSE) return "Synapse";
  else return "UNKNOWN-TYPE COMPARTMENT";
}
