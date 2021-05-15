import numpy as np
import matplotlib.pyplot as plt

alldata=np.genfromtxt("ode.txt")
plt.hist(alldata,200,histtype='step')
plt.xlabel('semi-major axis')
plt.ylabel('# of asteroids')
plt.show()