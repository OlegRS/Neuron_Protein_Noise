#include "Som_den_junction.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Soma.hpp"

Junction& Neuron::Som_den_junction::set_hopping_rate_constants() {

  fwd_mRNA_hop_rate = p_from->mRNA_diffusion_constant*3600/((p_from->length)*(p_from->length)) + p_from->mRNA_forward_trafficking_velocity*3600/p_from->length;
  bkwd_mRNA_hop_rate = p_to->mRNA_diffusion_constant*3600/((p_to->length)*(p_to->length)) + p_to->mRNA_backward_trafficking_velocity*3600/p_to->length;
  
  fwd_prot_hop_rate = p_from->protein_diffusion_constant*3600/((p_from->length)*(p_from->length)) + p_from->protein_forward_trafficking_velocity*3600/p_from->length;
  bkwd_prot_hop_rate = p_to->protein_diffusion_constant*3600/((p_to->length)*(p_to->length)) + p_to->protein_backward_trafficking_velocity*3600/p_to->length;

  return *this;
}
