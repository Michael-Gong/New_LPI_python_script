import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec

import multiprocessing as mp


######## Constant defined here ########
pi        =     3.1415926535897932384626
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     0.8e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 25,  
        }  
font2 = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 14,  
        }  
space_1 = 1
space_2 = 1
font_size = 25
marker_size=1.5
line_width=2.5
ave_span=200

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
    lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


if __name__ == '__main__':
  from_path = './PW_w020_trace/'
  to_path   = './PW_w020_trace_fig/'
  
  color_list= ['black','purple','royalblue','seagreen','darkorange','crimson'] 
  n_list    = [8,10,12,14,16]
  host=plt.subplot(1,1,1)
  plt.plot(np.zeros(201)+0.0,np.linspace(-50,50,201),':',color='k',linewidth=line_width)
  plt.plot(np.zeros(201)+0.3,np.linspace(-50,50,201),':',color='k',linewidth=line_width)
  for i in range(4):
      n=n_list[i]
      data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time1=header['time']
      ex = data['Electric Field/'+str.capitalize('ex')].data/exunit
      x  = data['Grid/Grid_mid'].data[0]/1.0e-6
#  host=plt.subplot(3,3,7)
      plt.plot(x,np.sum(ex[:,1500-ave_span:1500+ave_span],axis=1)/ave_span/2.,'-',color=color_list[i], label="t="+str(round(time1/1e-15,1))+' fs',linewidth=line_width,zorder=10-i)
  plt.legend(loc='best',fontsize=font_size,framealpha=0.5)
  plt.xlim(-0.3,1.1)
  plt.ylim(-5,15)
  plt.xlabel('$x\ [\mu m]$',fontdict=font)
  plt.ylabel('$E_x\ [m_ec\omega_0/|e|]$',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks([-5,0,5,10,15],fontsize=font_size);

  plt.subplots_adjust(left=0.15, bottom=0.17, right=0.92, top=0.97, wspace=0.01, hspace=0.01)
#      plt.title('Au51 at '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
  fig = plt.gcf()
  fig.set_size_inches(10, 5)
  fig.savefig('wrap_ex_time.png',format='png',dpi=160)
  plt.close("all")
  print('wrap_ex_time.png')
