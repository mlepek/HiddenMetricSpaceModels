
import numpy as np
import matplotlib.pyplot as plt


histogram = np.loadtxt("histo.txt", dtype=int)

x_axis = range(1, 100)

plt.loglog(x_axis, histogram[1:100])


plt.show()



