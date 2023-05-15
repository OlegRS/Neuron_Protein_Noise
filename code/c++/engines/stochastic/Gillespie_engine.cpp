#include <math.h>
#include "Gillespie_engine.hpp"

Compartment* Gillespie_engine::initialise_soma() {

  auto& soma = *p_neuron->p_soma;

  p_events[0] = &soma.gene_activation.set_rate(soma.gene_activation_rate*(soma.number_of_gene_copies-soma.n_active_genes));
  p_neuron->total_rate = soma.gene_activation.rate;
  p_events[1] = &soma.gene_deactivation.set_rate(soma.gene_deactivation_rate*soma.n_active_genes);
  p_neuron->total_rate += soma.gene_deactivation.rate;
  p_events[2] = &soma.mRNA_creation.set_rate(soma.transcription_rate*soma.n_active_genes);
  p_neuron->total_rate += soma.mRNA_creation.rate;
  p_events[3] = &soma.mRNA_decay.set_rate(soma.mRNA_decay_rate*soma.n_mRNAs);
  p_neuron->total_rate += soma.mRNA_decay.rate;
  p_events[4] = &soma.protein_creation.set_rate(soma.translation_rate*soma.n_mRNAs);
  p_neuron->total_rate += soma.protein_creation.rate*soma.n_mRNAs;
  p_events[5] = &soma.protein_decay.set_rate(soma.protein_decay_rate*soma.n_proteins);
  p_neuron->total_rate += soma.protein_decay.rate;

  ev_ind = 6;

  return &soma;
}

void Gillespie_engine::initialise_from(Compartment& comp) {

  for(auto& it_p_junc : comp.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    auto& par = *(p_junc->p_from); // Parent
    auto& desc = *(p_junc->p_to);  // Descendant

    if(p_junc->type() == DEN_SYN) {
      p_events[ev_ind++] = &p_junc->prot_hop_forward.set_rate(p_junc->fwd_prot_hop_rate*par.n_proteins);
      p_neuron->total_rate += p_junc->prot_hop_forward.rate;
      p_events[ev_ind++] = &p_junc->prot_hop_backward.set_rate(p_junc->bkwd_prot_hop_rate*desc.n_proteins);
      p_neuron->total_rate += p_junc->prot_hop_backward.rate;

      desc.mRNA_creation.rate = 0;
      desc.mRNA_decay.rate = 0;
      desc.protein_creation.rate = 0;
      desc.protein_decay.rate = 0;
    }
    else if(p_junc->type() == DEN_DEN || p_junc->type() == SOM_DEN) {
      p_events[ev_ind++] = &p_junc->mRNA_hop_forward.set_rate(p_junc->fwd_mRNA_hop_rate*par.n_mRNAs);
      p_neuron->total_rate += p_junc->mRNA_hop_forward.rate;
      p_events[ev_ind++] = &p_junc->mRNA_hop_backward.set_rate(p_junc->bkwd_mRNA_hop_rate*desc.n_mRNAs);
      p_neuron->total_rate += p_junc->mRNA_hop_backward.rate;
      p_events[ev_ind++] = &p_junc->prot_hop_forward.set_rate(p_junc->fwd_prot_hop_rate*par.n_proteins);
      p_neuron->total_rate += p_junc->prot_hop_forward.rate;
      p_events[ev_ind++] = &p_junc->prot_hop_backward.set_rate(p_junc->bkwd_prot_hop_rate*desc.n_proteins);
      p_neuron->total_rate += p_junc->prot_hop_backward.rate;

      desc.mRNA_creation.rate = 0;
      p_events[ev_ind++] = &desc.mRNA_decay.set_rate(desc.mRNA_decay_rate*desc.n_mRNAs);
      p_neuron->total_rate += desc.mRNA_decay.rate;
      p_events[ev_ind++] = &desc.protein_creation.set_rate(desc.translation_rate*desc.n_mRNAs);
      p_neuron->total_rate += desc.protein_creation.rate;
      p_events[ev_ind++] = &desc.protein_decay.set_rate(desc.protein_decay_rate*desc.n_proteins);
      p_neuron->total_rate += desc.protein_decay.rate;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
    
    initialise_from(*((*it_p_junc)->p_to));
  }
}

inline double Gillespie_engine::draw_delta_t() {
  double s = 0;
  for(auto& p_ev : p_events)
    s += p_ev->rate;

  // if(abs(p_neuron->total_rate-s)>.1) std::cerr << "abs(p_neuron->total_rate-s)>1\n";
  // std::cerr << "total_rate = " << p_neuron->total_rate << ", s = " << s << ", r = " << r << std::endl;

  p_neuron->total_rate = s; //Total rate should be computed without the above sum! Numerical stability is a problem
  
  // Sampling event time
  double delta_t = -1/p_neuron->total_rate*log(1-rnd());

  return delta_t; 
}

void Gillespie_engine::update_Gillespie() {
  // Sampling the event
  double r = rnd()*p_neuron->total_rate;

  size_t i=0;
  for(double sum=p_events[0]->rate; sum<r; sum += p_events[i]->rate)
    ++i;
  
  (*p_events[i])(); // Triggering the event
}

void Gillespie_engine::run_Gillespie(const double& time) {
  initialise_soma();
  initialise_from(*p_neuron->p_soma);

  // FOR PRELIMINARY TESTING ONLY!!!
  std::cout << "t," << "Soma_AG," << "Soma_mRNA,"<< "Soma_Prot,";
  for(auto& p_ds : p_neuron->p_dend_segments)
    std::cout << p_ds->name + "_mRNA," << p_ds->name+"_Prot,";
  for(auto& p_ds : p_neuron->p_synapses)
    std::cout << p_ds->name + "_Prot,";
  std::cout << std::endl;

  double t_prev=-1.1, delta_t_write = 1;
  for(double t = 0; t<time; t+=draw_delta_t()) {
    if(t-t_prev > delta_t_write) {
      std::cout << t << ',' << p_neuron->p_soma->n_active_genes << ',' << p_neuron->p_soma->n_mRNAs << ',' << p_neuron->p_soma->n_proteins << ',';
      for(auto& p_ds : p_neuron->p_dend_segments)
        std::cout << p_ds->n_mRNAs << ',' << p_ds->n_proteins << ',';
      for(auto& p_ds : p_neuron->p_synapses)
        std::cout << p_ds->n_proteins << ',';
      std::cout << std::endl;
      t_prev = t;
    }
    update_Gillespie();
  }
}

void Gillespie_engine::run_Gillespie(const std::list<double>& times) {
  initialise_soma();
  initialise_from(*p_neuron->p_soma);
  //for(double t = 0; t<time; t+=update_Gillespie());

  // FOR PRELIMINARY TESTING ONLY!!!
  std::cout << "t," << "time," << "Soma_AG," << "Soma_mRNA,"<< "Soma_Prot";
  for(auto& p_ds : p_neuron->p_dend_segments)
    std::cout << ',' + p_ds->name + "_mRNA," + p_ds->name + "_Prot";
  for(auto& p_ds : p_neuron->p_synapses)
    std::cout << ',' + p_ds->name + "_Prot";
  std::cout << std::endl;

  double t = 0;

  auto it_times=times.begin();
  while(t<=*it_times && it_times!=times.end()) { std::cerr << "A\n";
    std::cout << t << ',' << *it_times << ',' << p_neuron->p_soma->n_active_genes << ',' << p_neuron->p_soma->n_mRNAs << ',' << p_neuron->p_soma->n_proteins;
    for(auto& p_ds : p_neuron->p_dend_segments)
      std::cout << ',' << p_ds->n_mRNAs << ',' << p_ds->n_proteins;
    for(auto& p_ds : p_neuron->p_synapses)
      std::cout << ',' << p_ds->n_proteins;
    std::cout << std::endl;
    
    t+=draw_delta_t();
    while(t>*it_times && it_times!=times.end()) { std::cerr << "B\n";
      std::cout << t << ',' << *it_times << ',' << p_neuron->p_soma->n_active_genes << ',' << p_neuron->p_soma->n_mRNAs << ',' << p_neuron->p_soma->n_proteins;
      for(auto& p_ds : p_neuron->p_dend_segments)
        std::cout << ',' << p_ds->n_mRNAs << ',' << p_ds->n_proteins;
      for(auto& p_ds : p_neuron->p_synapses)
        std::cout << ',' << p_ds->n_proteins;
      std::cout << std::endl;
      ++it_times;
    }
    update_Gillespie();
    ++it_times;
  }



  // for(auto it_times=times.begin(); it_times!=times.end(); ++it_times) {
  //   std::cerr << "A\n";
  //   std::cout << t << ',' << *it_times << ',' << p_neuron->p_soma->n_active_genes << ',' << p_neuron->p_soma->n_mRNAs << ',' << p_neuron->p_soma->n_proteins;
  //   for(auto& p_ds : p_neuron->p_dend_segments)
  //     std::cout << ',' << p_ds->n_mRNAs << ',' << p_ds->n_proteins;
  //   for(auto& p_ds : p_neuron->p_synapses)
  //     std::cout << ',' << p_ds->n_proteins;
  //   std::cout << std::endl;

  //   while(t<=*it_times) {
  //     t+=draw_delta_t();
  //     while(t>*it_times && it_times!=times.end()) { std::cerr << "B\n";
  //       std::cout << t << ',' << *it_times << ',' << p_neuron->p_soma->n_active_genes << ',' << p_neuron->p_soma->n_mRNAs << ',' << p_neuron->p_soma->n_proteins;
  //       for(auto& p_ds : p_neuron->p_dend_segments)
  //         std::cout << ',' << p_ds->n_mRNAs << ',' << p_ds->n_proteins;
  //       for(auto& p_ds : p_neuron->p_synapses)
  //         std::cout << ',' << p_ds->n_proteins;
  //       std::cout << std::endl;
  //       ++it_times;
  //     }
  //     update_Gillespie();
  //   }
  // }
  
  // for(auto& time : times) {
  //   std::cout << t << ',' << time << ',' << p_neuron->p_soma->n_active_genes << ',' << p_neuron->p_soma->n_mRNAs << ',' << p_neuron->p_soma->n_proteins;
  //   for(auto& p_ds : p_neuron->p_dend_segments)
  //     std::cout << ',' << p_ds->n_mRNAs << ',' << p_ds->n_proteins;
  //   for(auto& p_ds : p_neuron->p_synapses)
  //     std::cout << ',' << p_ds->n_proteins;
  //   std::cout << std::endl;
    
  //   while(t<=time) {
  //     t+=draw_delta_t();
  //     update_Gillespie();
  //   }
  // }
}
