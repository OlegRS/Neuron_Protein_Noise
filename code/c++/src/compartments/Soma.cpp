#include "../../include/compartments/Soma.hpp"

std::ostream& operator<<(std::ostream& os, const Soma& som) {
  os << som.get_name() << ", out_junctions:\n";
  for(auto& it_p_out_junc : som.it_p_out_junctions)
    os << **it_p_out_junc << std::endl;
    
  return os;
}
