import numpy as np
import matplotlib.pyplot as plt

alldata=np.genfromtxt("odempi.txt")
plt.hist(alldata,200,histtype='step')
plt.xlabel('semi-major axis')
plt.ylabel('# of asteroids')
plt.title('histgram of the numbers of astroids with distance to Sun')
plt.show()
