import compartment.py
import soma.py
import dendrite.py
import dendritic_branch.py
import synapse.py

class neuron(compartment):
    # General parameters of the neuron
    tau1 = 1/(lambda1p+lambda1m)
    tau2 
    tau3

    # Constructors
    def __init__(self):
        
