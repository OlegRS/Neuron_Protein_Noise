#ifndef __DEN_DEN_JUNCTION_HPP__
#define __DEN_DEN_JUNCTION_HPP__

#include "../Neuron.hpp"

class Neuron::Den_den_junction : public Junction {
  friend class Analytic_engine;

  // // Parameters
  // double mRNA_diffusion_rate = .24*3600/(length*length); //  /hour; mRNA diffusion rate constant (3.4e-3 um2/s for CaMKII)
  // double mRNA_forward_trafficking_rate = .1*5e-2*3600/length; // /hour; mRNA forward trafficking rate  (4e-2 um/s for CaMKII net rate)
  // double mRNA_backward_trafficking_rate = .1*1e-2*3600/length; // /hour; mRNA forward trafficking rate  (4e-2 um/s for CaMKII net rate)
  // double protein_diffusion_rate = .24*3600/(length*length);
  // double protein_forward_trafficking_rate = 0;
  // double protein_backward_trafficking_rate = 0;
  
  // double protein_forward_hopping_rate;
  // double protein_backward_hopping_rate;
  // double mRNA_forward_hopping_rate;
  // double mRNA_backward_hopping_rate;
  
public:
  Den_den_junction();
  Den_den_junction(Compartment* p_from, Compartment* p_to) : Junction(p_from, p_to)
  { set_hopping_rate_constants(); }
  
  Junction& set_hopping_rate_constants();

  // std::vector<double> hopping_rates() {return std::vector<double>{fwd_mRNA_hop_rate, bkwd_mRNA_hop_rate, fwd_prot_hop_rate, bkwd_prot_hop_rate};}

  Junction::Type type() const {return DEN_DEN;}
};


#endif
