
import numpy as np
import matplotlib.pyplot as plt


d = np.loadtxt("out.txt", dtype=float)


hist, bin_edges = np.histogram(d, 1000)

x_axis = range(1, 100)

plt.plot((bin_edges[0:len(bin_edges)-1]+bin_edges[1:len(bin_edges)])/2, hist)


plt.show()



