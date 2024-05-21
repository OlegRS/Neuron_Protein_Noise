#include "Junction.hpp"

std::ostream& operator<<(std::ostream &os, const Junction &junc) {
  os << junc.p_from->name << "--[" << junc.type() << "]-->" << junc.p_to->name
     << " fwd_mRNA_hop_rate=" << junc.fwd_mRNA_hop_rate
     << " bkwd_mRNA_hop_rate=" << junc.bkwd_mRNA_hop_rate
     << " fwd_prot_hop_rate=" << junc.fwd_prot_hop_rate
     << " bkwd_prot_hop_rate=" << junc.bkwd_prot_hop_rate;
  return os;
}
