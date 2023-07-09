#ifndef __SOMA_HPP__
#define __SOMA_HPP__

#include "../Neuron.hpp"

class Soma : public Compartment {
  friend class Neuron::Som_den_junction;
  friend class Analytic_engine;
  friend class Gillespie_engine;
  
  // Parameters
  unsigned int number_of_gene_copies = 3; //For CaMKIIa it is 2 (Fonkeu)
  double gene_activation_rate = 1;
  double gene_deactivation_rate = 0;

  // For Monte Carlo engines
  struct Gene_activation : public Event {
    Gene_activation(Soma* p_loc) : Event(p_loc) {}
    Event::Type type() {return GENE_ACTIVATION;}
    void operator()();
  } gene_activation;

  struct Gene_deactivation : public Event {
    Gene_deactivation(Soma* p_loc) : Event(p_loc) {}
    Event::Type type() {return GENE_DEACTIVATION;}
    void operator()();
  } gene_deactivation;

  double n_active_genes_expectation = gene_activation_rate/(gene_activation_rate + gene_deactivation_rate)*number_of_gene_copies;
  unsigned int n_active_genes=0;
  
public:
  Soma(const std::string& name="no_name", const double& length=200) : Compartment(length, name), gene_activation(this), gene_deactivation(this) {
    transcription_rate = (3.*200/*dend_length*//10000)*0.001*3600; // /hour; mRNA transcription rate (0.001/s CaMKII Fonkeu) // THE FACTOR IN () ACCOUNTS FOR THE REDUCED LENGTH OF THE SIMPLE MODEL DENDRITE COMPARED TO THE REAL NEURONS
    mRNA_decay_rate = 1.2e-5*3600;
    translation_rate = 0.021*3600;
    protein_decay_rate = 1.21e-6*3600;
  }

  Soma(unsigned int &number_of_gene_copies, double &gene_activation_rate, double &gene_deactivation_rate, double &transcription_rate, double &mRNA_decay_rate, double &translation_rate, double &protein_decay_rate, unsigned int active_genes_number = 0, unsigned int protein_number = 0, unsigned int mRNA_number = 0);

  Soma& set_gene_activation_rate(const double& rate) {gene_activation_rate=rate; return *this;}
  Soma& set_gene_deactivation_rate(const double& rate) {gene_deactivation_rate = rate; return *this;}
  Soma& set_number_of_gene_copies(const unsigned int& N) {number_of_gene_copies = N; return *this;}
  Soma& set_transcription_rate(const double& rate) {transcription_rate=rate; return *this;};
  
  Compartment::Type type() const {return SOMA;}
  
  Compartment* operator+=(Compartment&);

  friend std::ostream& operator<<(std::ostream&, const Soma&);
 
  Compartment* add(Compartment&);
};

#endif
