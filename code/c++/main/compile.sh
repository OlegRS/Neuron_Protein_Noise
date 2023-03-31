#!/bin/bash

rm ../bin/exe/*

[ -d "../../data" ] || mkdir "../../data"

[ -d "../../bin" ] || mkdir "../../bin"

[ -d "../../bin/exe" ] || mkdir "../../bin/exe"


g++ ./main.cpp ../Neuron.cpp ../engines/analytic/Analytic_engine.cpp ../engines/stochastic/Gillespie_engine.cpp ../junctions/Junction.cpp ../junctions/Type.cpp ../junctions/Den_den_junction.cpp ../junctions/Som_den_junction.cpp ../junctions/Den_syn_junction.cpp ../compartments/Type.cpp ../compartments/Soma.cpp ../compartments/Synapse.cpp ../compartments/Dendritic_segment.cpp ../compartments/Compartment.cpp ../engines/stochastic/events/Event.cpp ../engines/stochastic/events/Protein_creation.cpp ../engines/stochastic/events/Protein_decay.cpp ../engines/stochastic/events/MRNA_creation.cpp ../engines/stochastic/events/MRNA_decay.cpp ../engines/stochastic/events/Gene_activation.cpp ../engines/stochastic/events/Gene_deactivation.cpp ../engines/stochastic/events/MRNA_hop_forward.cpp ../engines/stochastic/events/MRNA_hop_backward.cpp ../engines/stochastic/events/Prot_hop_forward.cpp ../engines/stochastic/events/Prot_hop_backward.cpp -o ../../bin/exe/main  -std=c++14 -O2 -larmadillo