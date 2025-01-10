import sys
sys.path.append('../build')

import SGEN_Py as sg

import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

Dendrie_length = 5000 #um
N_dendritic_segments = 1000

soma = sg.Soma("soma", Dendrie_length/N_dendritic_segments)

dendritic_segments = [sg.Dendritic_segment(soma, "d_1-1", Dendrie_length/N_dendritic_segments)]

for i in np.arange(1,N_dendritic_segments):
    dendritic_segments.append(sg.Dendritic_segment(dendritic_segments[i-1], "d_1-" + str(i+1), Dendrie_length/N_dendritic_segments))

neuron = sg.Neuron(soma, "Test_neuron")

ae = sg.Analytic_engine(neuron)

print("Computing mRNA expectations...")
mRNA_expectations = np.array(ae.stationary_mRNA_expectations())
print("Computing protein expectations...")
prot_expectations = np.array(ae.stationary_protein_expectations())
# print("Computing gene-mRNA correlations...")
# gene_mRNA_covariances = np.array(ae.stationary_gene_mRNA_covariances())
# print("Computing mRNA-mRNA correlations...")
# mRNA_mRNA_covariances = np.array(ae.stationary_mRNA_mRNA_covariances())
# print("Computing gene-protein correlations...")
# gene_prot_covariances = np.array(ae.stationary_gene_protein_covariances())
# print("Computing mRNA-protein correlations...")
# mRNA_prot_covariances = np.array(ae.stationary_mRNA_protein_covariances())
# print("Computing protein-protein correlations...")
# prot_prot_covariances = np.array(ae.stationary_protein_protein_covariances())

# Get the neuron segments and their volumes
me = sg.Morphologic_engine(neuron)
segments = me.segments()
volumes = me.volumes()

start_points = [segments[i][0][:3] for i in range(len(segments))]
end_points = [segments[i][1][:3] for i in range(len(segments))]
radii = [segments[i][0][3] for i in range(len(segments))]

prot_concentrations = prot_expectations/volumes
segment_values = np.flip(prot_concentrations)

# Convert the neuron morphology into a mesh for visualization in PyVista

# Make a tube for each segment to represent dendrites
tubes = []
all_scalars = [] # Collect scalars for the final mesh
for i in range(0,len(segments)):
    start, end = start_points[i], end_points[i]
    radius = radii[i]  # Radius of the current segment
    
    # Create a tube (cylinder) between two coordinates
    tube = pv.Line(start, end).tube(radius=radius)
    tubes.append(tube)

    # Extend scalars to match the number of points in the current tube
    all_scalars.append(np.full(tube.n_cells, segment_values[i]))

# Combine all tubes into a single mesh
neuron_mesh = tubes[0] if tubes else None
for tube in tubes[1:]:
    neuron_mesh += tube

# Concatenate scalars for all segments
flat_scalars = np.concatenate(all_scalars)

# Assign scalars to the combined mesh
neuron_mesh.cell_data["Protein Levels"] = flat_scalars

# Visualize the neuron
plotter = pv.Plotter()
plotter.add_mesh(neuron_mesh, scalars="Protein Levels", cmap="coolwarm", show_edges=False)
plotter.show_axes()
plotter.show()


# print(neuron)
