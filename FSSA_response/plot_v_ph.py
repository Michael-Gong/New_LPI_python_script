import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors
import scipy.ndimage as ndimage



######## Constant defined here ########
pi        =     3.1415926535897932384626
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     1.0e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 25,  
       }  

font2 = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 25,  
       }

font_size =25

to_path   = './two/'
grid_data    = np.load(to_path+'ey_vph.npy')
grid_time = np.load(to_path+'grid_time.npy')/1e-15
grid_x    = np.load(to_path+'grid_x.npy')

xx        =  np.load(to_path+'xx2d.npy')
px  =  np.load(to_path+'px2d.npy')
py  =  np.load(to_path+'py2d.npy')
gg  =  (1+px**2+py**2)**0.5
plt.subplot(1,1,1) 
#axin1 = inset_axes(ax, width='15%', height='5%', loc='upper left')
#axin2 = inset_axes(ax, width='15%', height='5%', loc='upper center')

X, Y = np.meshgrid(grid_x, grid_time) 
levels = np.linspace(-60.1, 60.1, 50)
grid_data[grid_data < -60]=-60
grid_data[grid_data >  60]= 60
plt.contourf(30+X-v0*Y*1e-15/1e-6, Y, grid_data.T, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.bwr)
cbar=plt.colorbar(ticks=np.linspace(-60,60,3))#,orientation="horizontal")
cbar.set_label('$E_y$ [$m_ec\omega/e$]', fontdict=font2)
for n in range(22):
    plt.plot(xx[n,:]+30-v0*grid_time*1e-15/1e-6, grid_time,'-k')
#plt.contourf(X, Y, grid_data.T, levels=levels, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()), cmap=cm.bwr)
n = 3
plt.scatter(xx[n,:]+30-v0*grid_time*1e-15/1e-6, grid_time, c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1500), s=10, cmap='rainbow', edgecolors='None', alpha=1,zorder=3)


#norm_x = matplotlib.colors.Normalize(vmin=np.min(workx2d[index,0:1299:20]),vmax=np.max(workx2d[index,0:1299:20]))
#print(np.shape(grid_t),np.shape(xx[index,0:1299:10]-(n-15)))
#plt.scatter(grid_t[5::2], xx[index,0:1299:20]-(grid_t[5::2]/3.3333333-15)-12, c=workx2d[index,0:1299:20], norm=norm_x, s=60, cmap='hot', edgecolors='black')
#cbar=plt.colorbar()#orientation="horizontal")
#cbar.set_label('Work$_x$ [$m_ec^2$]', fontdict=font2)

#plt.xlabel('Energy [MeV]',fontdict=font)
plt.ylabel('time [fs]',fontdict=font)
plt.xlabel('$\Delta=x-ct+30$ [$\mu$m]',fontdict=font)
plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#plt.yscale('log')
plt.ylim(0,320)
plt.xlim(-5,20)
#plt.legend(loc='best',fontsize=20,framealpha=0.5)
plt.grid()

fig = plt.gcf()
fig.set_size_inches(13.2, 10.3)
fig.savefig(to_path+'move_window_1.png',format='png',dpi=160)
plt.close("all")

