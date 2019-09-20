import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Data generation
alpha = np.linspace(1, 8, 5)
t = np.linspace(0, 5, 16)
T, A = np.meshgrid(t, alpha)
data = np.exp(-T * (1. / A))

# Plotting
fig = plt.figure()
ax = fig.gca(projection = '3d')

Xi = T.flatten()
Yi = A.flatten()
Zi = np.zeros(data.size)

dx = .25 * np.ones(data.size)
dy = .25 * np.ones(data.size)
dz = data.flatten()

ax.set_xlabel('T')
ax.set_ylabel('Alpha')
ax.bar3d(Xi, Yi, Zi, dx, dy, dz, color = 'w')

#plt.show()

plt.subplots_adjust(top=0.93, bottom=0.16, left=0.1, right=0.95, hspace=0.10, wspace=0.05)

fig = plt.gcf()
fig.set_size_inches(8, 7.2)
#fig.set_size_inches(5, 4.5)
fig.savefig('./3d_bar_plot.png',format='png',dpi=160)
plt.close("all")

  
          
