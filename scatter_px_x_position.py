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
wavelength=     1.0e-6
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
marker_size=0.5

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
  part_name = 'carbon'
  print('plotting '+part_name)
  from_path = './PW_w020/'
  to_path   = './PW_w020_fig/'
  px_d = np.loadtxt(from_path+part_name+'_px.txt')
  py_d = np.loadtxt(from_path+part_name+'_py.txt')
  xx_d = np.loadtxt(from_path+part_name+'_xx.txt')
  yy_d = np.loadtxt(from_path+part_name+'_yy.txt')
  start   =  1#14  # start time
  stop    =  19#28  # end time
  step    =  1  # the interval or step
  for n in range(start,stop+step,step):   
      data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time1=header['time']

      plt.subplot(2,1,1)
      plt.scatter(xx_d[:,n-start], px_d[:,n-start], c=xx_d[:,0], cmap='rainbow', norm=colors.Normalize(vmin=0, vmax=0.3) , s=marker_size, marker='.',alpha=0.8,label=part_name+' $p_x-x$',zorder=1,lw=0)
      cbar=plt.colorbar(pad=0.01,ticks=np.linspace(0, 0.3, 4))
      cbar.ax.tick_params(labelsize=font2['size']) 
      #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
      cbar.set_label('$x|_{t=0}$',fontdict=font2)        
      plt.text(4.,0.06,str(round(time1/1.0e-15,0))+' fs',fontdict=font2,color='k')
#      plt.legend(loc='best',fontsize=14)
      plt.xlim(-4,12)
      plt.ylim(-0.6,0.6)
      plt.xlabel('$x\ [\mu m]$',fontdict=font)
      plt.ylabel('$p_x\ [m_ic]$',fontdict=font)
      plt.xticks(fontsize=font_size); 
      plt.yticks(fontsize=font_size);

      plt.subplot(2,1,2)
      plt.scatter(xx_d[:,n-start], yy_d[:,n-start], c=xx_d[:,0], cmap='rainbow', norm=colors.Normalize(vmin=0, vmax=0.3) , s=marker_size, marker='.',alpha=0.8,label=part_name+' $x-y$',zorder=1,lw=0)
      cbar=plt.colorbar(pad=0.01,ticks=np.linspace(0, 0.3, 4))
      cbar.ax.tick_params(labelsize=font2['size']) 
      #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=font2['size'])        
      cbar.set_label('$x|_{t=0}$',fontdict=font2)        
      plt.text(4.,0.06,str(round(time1/1.0e-15,0))+' fs',fontdict=font2,color='k')
#      plt.legend(loc='best',fontsize=14)
      plt.xlim(-4,12)
      plt.ylim(-12,12)
      plt.xlabel('$x\ [\mu m]$',fontdict=font)
      plt.ylabel('$y\ [\mu m]$',fontdict=font)
      plt.xticks(fontsize=font_size); 
      plt.yticks(fontsize=font_size);

      plt.subplots_adjust(left=0.16, bottom=None, right=0.94, top=None, wspace=None, hspace=None)
#      plt.title('Au51 at '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
      fig = plt.gcf()
      fig.set_size_inches(12, 7.5)
      fig.savefig(to_path+part_name+'_trace_pos'+str(n).zfill(4)+'.png',format='png',dpi=320)
      plt.close("all")
      print('finised '+to_path+part_name+'_trace_pos'+str(n).zfill(4)+'.png')
