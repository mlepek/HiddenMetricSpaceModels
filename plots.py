



import numpy as np
import matplotlib.pyplot as plt




Ny = [300, 1000, 3000, 10000, 30000]

alfa = 5.0

gamma_30 = [0.505704, 0.428497, 0.345229, 0.288474, 0.225762]
gamma_28 = [0.547509, 0.494594, 0.47686, 0.456436, 0.426275]
gamma_25 = [0.545099, .537281, 0.544817, 0.556211, 0.542933]
gamma_22 = [0.33915, 0.356194, 0.367329, 0.3762, 0.371269]




plt.semilogx(Ny, gamma_22, 'rs-', Ny, gamma_25, 'g<-', Ny, gamma_28, 'm+-', Ny, gamma_30, 'c*-')

plt.title("alfa = 5.0")
plt.xlabel("wielkosc sieci, N")
plt.ylabel("success ratio, p_s")
plt.gca().legend(('gamma=2.2','gamma=2.5', 'gamma=2.8', 'gamma=3.0'))
plt.show()


N = 10000

gammy_do_alfa_15 = [2.2, 2.4, 2.6, 2.8]
gammy_do_alfa_50 = [2.2, 2.5, 2.8, 3.0]

alfa_15 = [0.266388, 0.416762, 0.475087, 0.499647]
alfa_50 = [0.3762, 0.556211, 0.456436, 0.288474]


plt.plot(gammy_do_alfa_15, alfa_15, 'rs-', gammy_do_alfa_50, alfa_50, 'g^-')
plt.show()