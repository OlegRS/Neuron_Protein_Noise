#include "../../include/junctions/Den_syn_junction.hpp"

Junction& Neuron::Den_syn_junction::set_hopping_rate_constants() {
  // fwd_prot_hop_rate = static_cast<Synapse*>(p_to)->protein_binding_rate;
  // bkwd_prot_hop_rate = static_cast<Synapse*>(p_to)->protein_unbinding_rate;

  double prot_fwd_diff_rate = p_from->protein_diffusion_constant*3600/((p_from->length)*(p_from->length));
  double prot_bkwd_diff_rate = p_to->protein_diffusion_constant*3600/((p_to->length)*(p_to->length));

  double prot_fwd_traff_rate = p_from->protein_forward_trafficking_velocity*3600/p_from->length;
  double prot_bkwd_traff_rate = p_to->protein_backward_trafficking_velocity*3600/p_to->length;

  fwd_prot_hop_rate = (prot_fwd_diff_rate + prot_fwd_traff_rate);// /static_cast<Dendritic_segment*>(p_from)->n_descending_DS;
  bkwd_prot_hop_rate = prot_bkwd_diff_rate + prot_bkwd_traff_rate;


  fwd_mRNA_hop_rate = 0;
  bkwd_mRNA_hop_rate = 0;

  return *this;
}
