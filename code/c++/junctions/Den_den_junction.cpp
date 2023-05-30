#include "Den_den_junction.hpp"
#include "../compartments/Dendritic_segment.hpp"

Junction& Neuron::Den_den_junction::set_hopping_rate_constants() {

  double mRNA_fwd_diff_rate = p_from->mRNA_diffusion_constant*3600/((p_from->length)*(p_from->length));
  double mRNA_bkwd_diff_rate = p_to->mRNA_diffusion_constant*3600/((p_to->length)*(p_to->length));

  double mRNA_fwd_traff_rate = p_from->mRNA_forward_trafficking_velocity*3600/p_from->length;
  double mRNA_bkwd_traff_rate = p_to->mRNA_backward_trafficking_velocity*3600/p_to->length;

  double prot_fwd_diff_rate = p_from->protein_diffusion_constant*3600/((p_from->length)*(p_from->length));
  double prot_bkwd_diff_rate = p_to->protein_diffusion_constant*3600/((p_to->length)*(p_to->length));

  double prot_fwd_traff_rate = p_from->protein_forward_trafficking_velocity*3600/p_from->length;
  double prot_bkwd_traff_rate = p_to->protein_backward_trafficking_velocity*3600/p_to->length;

    
  fwd_mRNA_hop_rate = (mRNA_fwd_diff_rate + mRNA_fwd_traff_rate)/static_cast<Dendritic_segment*>(p_from)->n_descending_DS;
  bkwd_mRNA_hop_rate = mRNA_bkwd_diff_rate + mRNA_bkwd_traff_rate;

  fwd_prot_hop_rate = (prot_fwd_diff_rate + prot_fwd_traff_rate)/static_cast<Dendritic_segment*>(p_from)->n_descending_DS;
  bkwd_prot_hop_rate = prot_bkwd_diff_rate + prot_bkwd_traff_rate;

  std::cerr << "p_from->length = " << p_from->length << ", p_to->length = " << p_to->length << std::endl;
  std::cerr << "DEN_DEN:\n"
            << "fwd_mRNA_hop_rate=" << fwd_mRNA_hop_rate << ", bkwd_mRNA_hop_rate=" << bkwd_mRNA_hop_rate << ", fwd_prot_hop_rate=" << fwd_prot_hop_rate << ", bkwd_prot_hop_rate=" << bkwd_prot_hop_rate << std::endl;

  return *this;
}
