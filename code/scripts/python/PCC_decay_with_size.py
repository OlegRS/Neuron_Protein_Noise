import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt


n_comps = [10, 50, 100, 150, 200]
pcc_soma_last_ds = [.785, .784, .766, .747, .728]
pcc_soma_first_ds = [.9885, .964, .931, .901, .872]

plt.plot(n_comps, pcc_soma_last_ds, marker='*')
plt.show()
plt.plot(n_comps, pcc_soma_first_ds, marker='*')
plt.show()
