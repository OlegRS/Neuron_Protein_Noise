import sys
sys.path.append('../build')

import SGEN_Py as sg

import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

Dendrie_length = 5000 #um
N_dendritic_segments = 10

soma = sg.Soma("soma", Dendrie_length/N_dendritic_segments)

dendritic_segments = [sg.Dendritic_segment(soma, "d_1-1", Dendrie_length/N_dendritic_segments)]

for i in np.arange(1,N_dendritic_segments):
    dendritic_segments.append(sg.Dendritic_segment(dendritic_segments[i-1], "d_1-" + str(i+1), Dendrie_length/N_dendritic_segments))

neuron = sg.Neuron(soma, "Test_neuron")

ae = sg.Analytic_engine(neuron, True)

ae.stationary_expectations_and_correlations()

# print(neuron)
