import sys
sys.path.append('../build')

import SGEN_Py as sg

import pyvista as pv
import numpy as np

Dendrie_length = 5000 #um
N_dendritic_segments = 10

soma = sg.Soma("soma", Dendrie_length/N_dendritic_segments, x=0, y=0, z=0, radius=50)

main_branch = [sg.Dendritic_segment(attached_to=soma,
                                    name = "d_1-1",
                                    length = Dendrie_length/N_dendritic_segments)]
for i in np.arange(1,N_dendritic_segments):
    main_branch.append(sg.Dendritic_segment(attached_to=main_branch[i-1],
                                            name="d_1-" + str(i+1),
                                            length=Dendrie_length/N_dendritic_segments))

secondary_branch = [sg.Dendritic_segment(attached_to=main_branch[5],#soma
                                         name="d_1_1-0",
                                         length=Dendrie_length/N_dendritic_segments,
                                         radius=5,
                                         d_theta=30*np.pi/360,
                                         d_phi=0)]
for i in np.arange(1,N_dendritic_segments):
    secondary_branch.append(sg.Dendritic_segment(attached_to=secondary_branch[i-1],
                                                 name="d_1_1-" + str(i+1),
                                                 length=Dendrie_length/N_dendritic_segments))

neuron = sg.Neuron(soma, "Test_neuron")

me = sg.Morphologic_engine(neuron)

segments = me.segments()

ae = sg.Analytic_engine(neuron)
ae.stationary_expectations_and_correlations()

# Extract the neuron segments and nodes
start_points = [segments[i][0][:3] for i in range(len(segments))]
end_points = [segments[i][1][:3] for i in range(len(segments))]
radii = [segments[i][0][3] for i in range(len(segments))]

# prot_expectations = np.genfromtxt("protein_expectations", delimiter='\n')
prot_expectations = np.genfromtxt("protein_expectations.dat", delimiter='\n')
segment_values = np.flip(np.log(prot_expectations))

# Convert the neuron morphology into a mesh for visualization in PyVista

# You can create a line or tube representation of the morphology by using segments (radii)
# Make a tube for each segment to represent dendrites
tubes = []
all_scalars = [] # Collect scalars for the final mesh
for i in range(len(segments)):
    start, end = start_points[i], end_points[i]
    radius = radii[i]  # Radius of the current segment
    print("segment_values[i]=", segment_values[i])
    print("radius=", radius)

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

    
# # Visualize the neuron
# plotter = pv.Plotter()
# plotter.add_mesh(neuron_mesh, scalars=all_scalars, show_edges=False, cmap="coolwarm")

# plotter.show_axes()

# plotter.show()
