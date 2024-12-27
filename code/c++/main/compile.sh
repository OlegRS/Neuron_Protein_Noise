#!/bin/bash

rm ../bin/exe/*

[ -d "../../data" ] || mkdir "../../data"

[ -d "../../bin" ] || mkdir "../../bin"

[ -d "../../bin/exe" ] || mkdir "../../bin/exe"


# g++ ./main.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../engines/stochastic/Gillespie_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/main  -std=c++14 -O2 -larmadillo

# g++ ./gillespie.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../engines/stochastic/Gillespie_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/gillespie  -std=c++14 -O2 -larmadillo

# g++ ./loading_from_swc.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../engines/stochastic/Gillespie_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/loading_from_swc  -std=c++14 -O2 -larmadillo

# g++ ./time_correlations.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/time_correlations  -std=c++14 -O2 -larmadillo

# g++ ./time_correlations_with_MC.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../engines/stochastic/Gillespie_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/time_correlations_with_MC  -std=c++14 -O2 -larmadillo

# g++ ./heterosynaptic_plasticity.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../engines/stochastic/Gillespie_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/heterosynaptic_plasticity  -std=c++14 -O2 -larmadillo

# g++ ./heterosynaptic_plasticity_MC.cpp ../Neuron.cpp ../engines/stochastic/Gillespie_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/heterosynaptic_plasticity_MC  -std=c++14 -O2 -larmadillo

# g++ ./susceptibilities.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp  ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/susceptibilities  -std=c++14 -O2 -larmadillo

# g++ ./full_neuron_HP_num.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp  ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/full_neuron_HP_num  -std=c++14 -O2 -larmadillo

# g++ ./heterosynaptic_plasticity_num.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp  ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/heterosynaptic_plasticity_num  -std=c++14 -O2 -larmadillo


# g++ ./heterosynaptic_plasticity_MC.cpp ../Neuron.cpp ../engines/stochastic/Gillespie_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/heterosynaptic_plasticity_MC  -std=c++14 -O2 -larmadillo


# g++ ./stationary_moments_for_cosyne.cpp ../Neuron.cpp ../engines/stochastic/Gillespie_engine.cpp ../engines/analytic/Analytic_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/stationary_moments_for_cosyne  -std=c++14 -O2 -larmadillo


# g++ ./susceptibilities.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp  ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/susceptibilities  -std=c++14 -O2 -larmadillo

# g++ ./stationary_moments_theor.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp  ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/stationary_moments_theor  -std=c++14 -O2 -larmadillo

# g++ ./properly_discretised_simple_neuron.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/properly_discretised_simple_neuron  -std=c++14 -O2 -larmadillo

# g++ ./properly_discretised_linear_neuron.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/properly_discretised_linear_neuron  -std=c++14 -O2 -larmadillo

# g++ ./pcc_densities_for_different_discretisations.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/pcc_densities_for_different_discretisations  -std=c++14 -O2 -larmadillo

# g++ ./semiinfinite_dendrite.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/semiinfinite_dendrite  -std=c++14 -O2 -larmadillo

g++ ./semiinfinite_dendrite_ODE_tests.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/semiinfinite_dendrite_ODE  -std=c++14 -O2 -larmadillo

# g++ ./semiinfinite_dendrite_ODE_sem_test.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/semiinfinite_dendrite_ODE_sem_test  -std=c++14 -O2 -larmadillo
