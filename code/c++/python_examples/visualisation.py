import sys
sys.path.append('../build')

import SGEN_Py as sg

import pyvista as pv
import neurom
import numpy as np

# Load the neuron morphology from the SWC file
# neuron = neurom.load_morphology("10_2REDO-850-GM18-Ctl-Ctl-Chow-BNL16A-CA1_Finished2e.swc")
neuron = neurom.load_morphology("DD13-10-c5-1.CNG_.swc")

# Extract the neuron segments and nodes
start_points = [neuron.segments[i][0][:3] for i in range(len(neuron.segments))]
end_points = [neuron.segments[i][1][:3] for i in range(len(neuron.segments))]
radii = [neuron.segments[i][1][3] for i in range(len(neuron.segments))]

# Convert the neuron morphology into a mesh for visualization in PyVista
# Create a mesh for the dendrites
mesh = pv.PolyData(coordinates)

segment_values = np.random.rand(len(neuron.segments))  # Example: random values

# You can create a line or tube representation of the morphology by using segments (radii)
# Make a tube for each segment to represent dendrites
dendrite_tubes = []
all_scalars = [] # Collect scalars for the final mesh
for i in range(len(neuron.segments)):
    start, end = start_points[i], end_points[i]
    radius = radii[i]  # Radius of the current segment

    # Create a tube (cylinder) between two coordinates
    tube = pv.Line(start, end).tube(radius=radius * 10)
    dendrite_tubes.append(tube)

    # Extend scalars to match the number of points in the current tube
    all_scalars.append(np.full(tube.n_cells, segment_values[i]))

# Combine all tubes into a single mesh
neuron_mesh = dendrite_tubes[0] if dendrite_tubes else None
for tube in dendrite_tubes[1:]:
    neuron_mesh += tube

# Visualize the neuron
plotter = pv.Plotter()
plotter.add_mesh(neuron_mesh, scalars=all_scalars, cmap="coolwarm", show_edges=False)

plotter.show_axes()

plotter.show()
