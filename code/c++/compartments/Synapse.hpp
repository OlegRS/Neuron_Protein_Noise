#ifndef __SYNAPSE_HPP__
#define __SYNAPSE_HPP__

#include "../Neuron.hpp"

class Synapse : public Compartment {
  friend class Neuron::Den_syn_junction;
  friend class Analytic_engine;
  
  //Parameters
  double protein_binding_rate = .6; // per hour
  double protein_unbinding_rate = 6;

  double length = 200/60,//2.5, //um
    protein_diffusion_constant = .24, //0.24 um2/s for CaMKII
    protein_forward_trafficking_velocity = 0, 
    protein_backward_trafficking_velocity = 0;
  
public:
  Synapse();
  Synapse(Compartment &parent, const std::string& name = "no_name");
  Synapse(Compartment &parent, const std::string& name, const double &protein_binding_rate, const double &protein_unbinding_rate, const double &protein_decay_rate=0, const unsigned int &protein_number = 0);

  Synapse& copy() const;

  const double& get_protein_binding_rate() const {return protein_binding_rate;}
  const double& get_protein_unbinding_rate() const {return protein_unbinding_rate;}

  Synapse& set_protein_binding_rate(const double& rate) {protein_binding_rate=rate; return *this;}
  Synapse& set_protein_unbinding_rate(const double& rate) {protein_unbinding_rate=rate; return *this;}

  Compartment::Type type() const {return SYNAPSE;}

  friend std::ostream& operator<<(std::ostream&, const Synapse&);
};

#endif
