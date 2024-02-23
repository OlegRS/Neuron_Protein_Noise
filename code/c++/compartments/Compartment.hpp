#ifndef __COMPARTMENT_HPP__
#define __COMPARTMENT_HPP__

#include <iostream>
#include <list>
#include "../engines/stochastic/events/Event.hpp"

class Neuron;
class Synapse;
class Dendritic_segment;
class Junction;
class Analytic_engine;
class Gillespie_engine;

class Compartment {// Abstract class
  friend class Neuron;
  friend class Synapse;
  friend class Dendritic_segment;
  friend class Analytic_engine;
  friend class Gillespie_engine;
  friend class Junction;
protected:
  
  struct Type {
#define SOMA 1 // Better replace these macros with enum
#define BASAL_DENDRITE 3
#define APICAL_DENDRITE 4
#define SYNAPSE 7

    short id; // Type id consistent with .swc neuron morphology format
  
    Type(const int& type_id) : id(type_id) {}
    operator std::string() const;

    friend bool operator==(const Type& type, const short&& id) {return type.id == id;}
    friend bool operator!=(const Type& type, const short&& id) {return type.id != id;}
    friend std::ostream& operator<<(std::ostream& os, const Type& type) {
      return os << std::string(type);
    }
  };

  std::string name;

  double transcription_rate = (3.*200/*dend_length*//10000) * .001*3600, // /hour; mRNA transcription rate (0.001/s CaMKII Fonkeu) // THE FACTOR IN () ACCOUNTS FOR THE REDUCED LENGTH OF THE SIMPLE MODEL DENDRITE COMPARED TO THE REAL NEURONS
    mRNA_decay_rate = 1.2e-5*3600,
    translation_rate = 0.021*3600,
    protein_decay_rate = 1.21e-6*3600,
    mRNA_diffusion_constant = 3.4e-3, // 3.4e-3 um2/s for CaMKII
    protein_diffusion_constant = .24, //0.24 um2/s for CaMKII
    mRNA_forward_trafficking_velocity = .5e-2, // 4e-2 um/s for CaMKII
    mRNA_backward_trafficking_velocity = .1e-2, // 4e-2 um/s for CaMKII
    protein_forward_trafficking_velocity = 0, 
    protein_backward_trafficking_velocity = 0;

  double x, y, z, r,
    length = 200; //micrometers

  std::list<Compartment*> p_descendants; // Descendants of the compartment in the tree

  std::list<std::list<Junction*>::iterator> it_p_out_junctions, // Outgoing junctions from compartment  
                                            it_p_in_junctions; // Incomming junctions to compartment
  
  Neuron *p_neuron = NULL; // Pointer to its neuron
  std::list<Compartment*>::iterator iterator; // Its position in its neuron's container

  // For Analytic_engine
  size_t o1_index, mRNA_ind, prot_ind, id;

  // For stochastic engines
  struct Protein_creation : public Event {
    Protein_creation(Compartment* p_loc) : Event(p_loc) {}
    Event::Type type() {return PROTEIN_CREATION;}
    void operator()();
  } protein_creation;
  struct Protein_decay : public Event {
    Protein_decay(Compartment* p_loc) : Event(p_loc) {}
    Event::Type type() {return PROTEIN_DECAY;}
    void operator()();
  } protein_decay;
  struct MRNA_creation : public Event {
    MRNA_creation(Compartment* p_loc) : Event(p_loc) {}
    Event::Type type() {return MRNA_CREATION;}
    void operator()();
  } mRNA_creation;
  struct MRNA_decay : public Event {
    MRNA_decay(Compartment* p_loc) : Event(p_loc) {}
    Event::Type type() {return MRNA_DECAY;}
    void operator()();
  } mRNA_decay;

  double n_mRNA_expectation, n_prot_expectation;
  unsigned int n_mRNAs=0, n_proteins=0;

  Compartment& clear_junctions() {
    it_p_out_junctions.clear();
    it_p_in_junctions.clear();
    return *this;
  }

public:
  Compartment(const std::string& name = "no_name") : name(name),
                                                     protein_creation(this),
                                                     protein_decay(this),
                                                     mRNA_creation(this),
                                                     mRNA_decay(this) {}
  
  Compartment(const double& length, const std::string& name = "no_name") : name(name),
                                                                           length(length),
                                                                           protein_creation(this),
                                                                           protein_decay(this),
                                                                           mRNA_creation(this),
                                                                           mRNA_decay(this) {}


  // Compartment(Compartment *parent, const std::string& name = "no_name") = 0; //This only needs to be implemented in actual compartments
  Compartment(const Compartment&);

  virtual Type type() const = 0;
  std::string get_name() const {return name;}

  Neuron* which_neuron() {return p_neuron;}

  Compartment& connect_to(Compartment&); // Linking another compartment
  Compartment& disconnect_from(Compartment&); // Unlinking another compartment

  Compartment& set_transcription_rate(const double& rate) {transcription_rate=rate; return *this;}
  Compartment& set_mRNA_decay_rate(const double& rate) {mRNA_decay_rate=rate; return *this;}
  Compartment& set_translation_rate(const double& rate) {translation_rate=rate; return *this;}
  Compartment& set_protein_decay_rate(const double& rate) {protein_decay_rate=rate; return *this;}

  ~Compartment() {};

  friend std::ostream& operator<<(std::ostream&, const Compartment&);
  friend std::ostream& operator<<(std::ostream&, const Junction&);
};

#endif
