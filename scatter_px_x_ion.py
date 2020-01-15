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
        'size'   : 8,  
        }  
space_1 = 1
space_2 = 1
font_size = 25
marker_size=0.02

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
  start   =  1#14  # start time
  stop    =  19#28  # end time
  step    =  1  # the interval or step

  color_list = ['blue','limegreen','red'] 

  for n in range(start,stop+step,step):   
      from_path = './PW_w020/'
      to_path   = './PW_w020_fig/'
  
      data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time1=header['time']
      fig,host = plt.subplots()
      if 'Grid/Particles/electron' in data.keys():
          ion_x = data['Grid/Particles/electron'].data[0]/1e-6
          ion_px = data['Particles/Px/electron'].data/(1000*m0*v0)
          plt.scatter(ion_x[::space_2], ion_px[::space_2], c=color_list[0], s=marker_size, marker='.',alpha=0.8,label='electron $p_x/1000m_ec$',zorder=1,lw=0)
      if 'Grid/Particles/carbon' in data.keys():
          ion_x = data['Grid/Particles/carbon'].data[0]/1e-6
          ion_px = data['Particles/Px/carbon'].data/(12.0*1836*m0*v0)
          plt.scatter(ion_x[::space_2], ion_px[::space_2], c=color_list[1], s=marker_size, marker='.',alpha=0.8,label='carbon $p_x/m_ic$',zorder=2,lw=0)
      if 'Grid/Particles/proton' in data.keys():
          ion_x = data['Grid/Particles/proton'].data[0]/1e-6
          ion_px = data['Particles/Px/proton'].data/(1836*m0*v0)
          plt.scatter(ion_x[::space_2], ion_px[::space_2], c=color_list[2], s=marker_size, marker='.',alpha=0.8,label='proton $p_x/m_ic$',zorder=3,lw=0)
      print('hehe1')
      plt.text(4.,0.06,str(round(time1/1.0e-15,0))+' fs',fontdict=font2,color='k')
      plt.legend(loc='best',fontsize=14)
      plt.xlim(-4,12)
      plt.ylim(-0.5,0.5)
      plt.xlabel('$x\ [\mu m]$',fontdict=font)
      plt.ylabel('$p_x\ [m_ic]$',fontdict=font)
      plt.xticks(fontsize=font_size); 
      plt.yticks(fontsize=font_size);
      plt.subplots_adjust(left=0.16, bottom=None, right=0.94, top=None, wspace=None, hspace=None)
#      plt.title('Au51 at '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
      par1 = host.twinx()
  #par2 = host.twinx()
  #par3 = host.twinx()
      ex = data['Electric Field/'+str.capitalize('ex')].data/exunit
      ey = data['Electric Field/'+str.capitalize('ey')].data/exunit
      x  = data['Grid/Grid_mid'].data[0]/1.0e-6
      p1, = par1.plot(x,np.sum(ex[:,1100:1300],axis=1)/200., '-',color='lime', label="$E_x$",linewidth=1)
      p1, = par1.plot(x,np.sum(ey[:,1100:1300],axis=1)/200., '-',color='k', label="$E_y$",linewidth=1)
#
      par1.legend(loc='upper center',fontsize=20,framealpha=0.0)
      par1.set_ylim(-40,40)
      par1.set_ylabel('$E_x$ [$m_ec\omega/|e|$]',fontdict=font,color='lime')
      par1.tick_params(axis='y',labelsize=25,colors='lime')
      fig = plt.gcf()
      fig.set_size_inches(12, 7.5)
      fig.savefig(to_path+'scatter_px_x_'+str(n).zfill(4)+'.png',format='png',dpi=160)
      plt.close("all")
      print('finised '+str(n).zfill(4))
