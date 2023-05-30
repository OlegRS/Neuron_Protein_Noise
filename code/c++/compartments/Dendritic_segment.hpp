#ifndef __DENDRITIC_SEGMENT_HPP__
#define __DENDRITIC_SEGMENT_HPP__

#include "../Neuron.hpp"

class Dendritic_segment : public Compartment {
  friend class Neuron::Den_den_junction;
  friend class Neuron::Som_den_junction;
  friend class Analytic_engine;

  size_t n_descending_DS = 0;// Number of decstnding dendritic segments 
  
  //Parameters
  
public:
  Dendritic_segment(const double &mRNA_diffusion_rate, const double &mRNA_forward_trafficking_rate, const double &protein_diffusion_rate, const double &protein_trafficking_rate, const double &translation_rate, const unsigned int &mRNA_numer = 0, const unsigned int &protein_number = 0);
  Dendritic_segment(Compartment &parent, const std::string& name = "no_name", const double& length=200);

  Dendritic_segment& copy() const;

  Compartment::Type type() const {return APICAL_DENDRITE;}

  friend std::ostream& operator<<(std::ostream&, const Dendritic_segment&);
};

#endif
