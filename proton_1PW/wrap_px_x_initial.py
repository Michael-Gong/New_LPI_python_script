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
ave_span=100

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
  part_name = 'proton'
  print('plotting '+part_name)
  from_path = './PW_w020/'
  to_path   = './PW_w020_fig/'
  px_d = np.loadtxt(from_path+part_name+'_px.txt')
  py_d = np.loadtxt(from_path+part_name+'_py.txt')
  xx_d = np.loadtxt(from_path+part_name+'_xx.txt')
  yy_d = np.loadtxt(from_path+part_name+'_yy.txt')
  start   =  1#14  # start time
  stop    =  16#28  # end time
  step    =  1  # the interval or step
      
  n=1
  data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time1=header['time']
#  ex = data['Electric Field/'+str.capitalize('ex')].data/exunit
#  den_e = data['Derived/Number_Density/electron'].data/denunit
#  den_p = data['Derived/Number_Density/proton'].data/denunit
#  den_c = data['Derived/Number_Density/carbon'].data/denunit
#  x  = data['Grid/Grid_mid'].data[0]/1.0e-6
#  host=plt.subplot(3,3,7)
#  plt.plot(x,np.sum(den_p[:,1500-ave_span:1500+ave_span],axis=1)/ave_span/2.+np.sum(den_c[:,1500-ave_span:1500+ave_span],axis=1)/ave_span/2.*6,'-',color='crimson', label=r"$\rho_i$",linewidth=line_width)
#  plt.plot(x,np.sum(den_e[:,1500-ave_span:1500+ave_span],axis=1)/ave_span/2.,'--',color='limegreen', label=r"$\rho_e$",linewidth=line_width)
#  plt.xlim(-0.5,4.5)
#  plt.ylim(-5,65)
#  plt.xlabel('$x\ [\mu m]$',fontdict=font)
#  plt.ylabel(r'$\rho_e\ or\ \rho_i\ [|e|n_c]$',fontdict=font)
#  plt.xticks(fontsize=font_size); 
#  plt.yticks(fontsize=font_size);
#  par1 = host.twinx()
#  p1, = par1.plot(x,np.sum(ex[:,1250:1750],axis=1)/500.,'-',color='blue', label="$E_x$",linewidth=line_width)
#  par1.set_ylim(-10,10)
##  par1.set_ylabel('$E_x$ [$m_ec\omega/|e|$]',fontdict=font,color='blue')
#  par1.tick_params(axis='y',labelsize=0.001,colors='blue')
  plt.subplot(2,3,4)
  plt.scatter(xx_d[:,n-start], px_d[:,n-start], c=xx_d[:,0], cmap='rainbow', norm=colors.Normalize(vmin=0, vmax=0.3) , s=marker_size, marker='.',alpha=0.8,lw=0)
#  cbar=plt.colorbar(pad=0.01,ticks=np.linspace(0, 0.3, 4))
#  cbar.ax.tick_params(labelsize=font2['size']) 
#  cbar.set_label('$x|_{t=0}$',fontdict=font2)        
  plt.text(2.,0.55,str(round(time1/1.0e-15,0))+' fs',fontdict=font2,color='k')
  plt.xlim(-0.5,4.5)
  plt.ylim(-0.15,0.65)
  plt.xlabel('$x\ [\mu m]$',fontdict=font)
  plt.ylabel('$p_x\ [m_ic]$',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks([0,0.2,0.4,0.6],fontsize=font_size);
  plt.subplot(2,3,1)
  plt.scatter(xx_d[:,n-start], yy_d[:,n-start], c=xx_d[:,0], cmap='rainbow', norm=colors.Normalize(vmin=0, vmax=0.3) , s=marker_size, marker='.',alpha=0.8,label=part_name+' $x-y$',zorder=1,lw=0)
#  cbar=plt.colorbar(pad=0.01,ticks=np.linspace(0, 0.3, 4))
#  cbar.ax.tick_params(labelsize=font2['size']) 
#  cbar.set_label('$x|_{t=0}$',fontdict=font2)        
  plt.xlim(-0.5,4.5)
  plt.ylim(-7,7)
#  plt.xlabel('$x\ [\mu m]$',fontdict=font)
  plt.ylabel('$y\ [\mu m]$',fontdict=font)
  plt.xticks(fontsize=0.001); 
  plt.yticks(fontsize=font_size);

  n=5
  data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time1=header['time']
#  ex = data['Electric Field/'+str.capitalize('ex')].data/exunit
#  den_e = data['Derived/Number_Density/electron'].data/denunit
#  den_p = data['Derived/Number_Density/proton'].data/denunit
#  den_c = data['Derived/Number_Density/carbon'].data/denunit
#  x  = data['Grid/Grid_mid'].data[0]/1.0e-6
#  host=plt.subplot(3,3,8)
#  plt.plot(x,np.sum(den_p[:,1500-ave_span:1500+ave_span],axis=1)/ave_span/2.+np.sum(den_c[:,1500-ave_span:1500+ave_span],axis=1)/ave_span/2.*6,'-',color='crimson', label=r"$\rho_i$",linewidth=line_width)
#  plt.plot(x,np.sum(den_e[:,1500-ave_span:1500+ave_span],axis=1)/ave_span/2.,'--',color='limegreen', label=r"$\rho_e$",linewidth=line_width)
#  plt.xlim(-0.5,4.5)
#  plt.ylim(-5,65)
#  plt.xlabel('$x\ [\mu m]$',fontdict=font)
##  plt.ylabel(r'$\rho_e\ or\ \rho_i\ [|e|n_c]$',fontdict=font)
#  plt.xticks(fontsize=font_size); 
#  plt.yticks(fontsize=0.001);
#  par1 = host.twinx()
#  p1, = par1.plot(x,np.sum(ex[:,1250:1750],axis=1)/500.,'-',color='blue', label="$E_x$",linewidth=line_width)
#  par1.set_ylim(-10,10)
##  par1.set_ylabel('$E_x$ [$m_ec\omega/|e|$]',fontdict=font,color='blue')
#  par1.tick_params(axis='y',labelsize=0.001,colors='blue')
  plt.subplot(2,3,5)
  plt.scatter(xx_d[:,n-start], px_d[:,n-start], c=xx_d[:,0], cmap='rainbow', norm=colors.Normalize(vmin=0, vmax=0.3) , s=marker_size, marker='.',alpha=0.8,lw=0)
#  cbar=plt.colorbar(pad=0.01,ticks=np.linspace(0, 0.3, 4))
#  cbar.ax.tick_params(labelsize=font2['size']) 
#  cbar.set_label('$x|_{t=0}$',fontdict=font2)        
  plt.text(2.,0.55,str(round(time1/1.0e-15,0))+' fs',fontdict=font2,color='k')
  plt.xlim(-0.5,4.5)
  plt.ylim(-0.15,0.65)
  plt.xlabel('$x\ [\mu m]$',fontdict=font)
#  plt.ylabel('$p_x\ [m_ic]$',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks([0,0.2,0.4,0.6],fontsize=0.001);
  plt.subplot(2,3,2)
  plt.scatter(xx_d[:,n-start], yy_d[:,n-start], c=xx_d[:,0], cmap='rainbow', norm=colors.Normalize(vmin=0, vmax=0.3) , s=marker_size, marker='.',alpha=0.8,label=part_name+' $x-y$',zorder=1,lw=0)
#  cbar=plt.colorbar(pad=0.01,ticks=np.linspace(0, 0.3, 4))
#  cbar.ax.tick_params(labelsize=font2['size']) 
#  cbar.set_label('$x|_{t=0}$',fontdict=font2)        
  plt.xlim(-0.5,4.5)
  plt.ylim(-7,7)
#  plt.xlabel('$x\ [\mu m]$',fontdict=font)
#  plt.ylabel('$y\ [\mu m]$',fontdict=font)
  plt.xticks(fontsize=0.001); 
  plt.yticks(fontsize=0.001);


  n=7
  data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time1=header['time']
#  ex = data['Electric Field/'+str.capitalize('ex')].data/exunit
#  den_e = data['Derived/Number_Density/electron'].data/denunit
#  den_p = data['Derived/Number_Density/proton'].data/denunit
#  den_c = data['Derived/Number_Density/carbon'].data/denunit
#  x  = data['Grid/Grid_mid'].data[0]/1.0e-6
#  host=plt.subplot(3,3,9)
#  plt.plot(x,np.sum(den_p[:,1500-ave_span:1500+ave_span],axis=1)/ave_span/2.+np.sum(den_c[:,1500-ave_span:1500+ave_span],axis=1)/ave_span/2.*6,'-',color='crimson', label=r"$\rho_i$",linewidth=line_width)
#  plt.plot(x,np.sum(den_e[:,1500-ave_span:1500+ave_span],axis=1)/ave_span/2.,'--',color='limegreen', label=r"$\rho_e$",linewidth=line_width)
#  plt.xlim(-0.5,4.5)
#  plt.ylim(-5,65)
#  plt.xlabel('$x\ [\mu m]$',fontdict=font)
##  plt.ylabel(r'$\rho_e\ or\ \rho_i\ [|e|n_c]$',fontdict=font)
#  plt.xticks(fontsize=font_size); 
#  plt.yticks(fontsize=0.001);
#  par1 = host.twinx()
#  p1, = par1.plot(x,np.sum(ex[:,1250:1750],axis=1)/500.,'-',color='blue', label="$E_x$",linewidth=line_width)
#  par1.set_ylim(-10,10)
#  par1.set_ylabel('$E_x$ [$m_ec\omega/|e|$]',fontdict=font,color='blue')
#  par1.tick_params(axis='y',labelsize=font_size,colors='blue')
  plt.subplot(2,3,6)
  plt.scatter(xx_d[:,n-start], px_d[:,n-start], c=xx_d[:,0], cmap='rainbow', norm=colors.Normalize(vmin=0, vmax=0.3) , s=marker_size, marker='.',alpha=0.8,lw=0)
#  cbar=plt.colorbar(pad=0.01,ticks=np.linspace(0, 0.3, 4))
#  cbar.ax.tick_params(labelsize=font2['size']) 
#  cbar.set_label('$x|_{t=0}$',fontdict=font2)        
  plt.text(2.,0.55,str(round(time1/1.0e-15,0))+' fs',fontdict=font2,color='k')
  plt.xlim(-0.5,4.5)
  plt.ylim(-0.15,0.65)
  plt.xlabel('$x\ [\mu m]$',fontdict=font)
#  plt.ylabel('$p_x\ [m_ic]$',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks([0,0.2,0.4,0.6],fontsize=0.001);
  plt.subplot(2,3,3)
  plt.scatter(xx_d[:,n-start], yy_d[:,n-start], c=xx_d[:,0], cmap='rainbow', norm=colors.Normalize(vmin=0, vmax=0.3) , s=marker_size, marker='.',alpha=0.8,label=part_name+' $x-y$',zorder=1,lw=0)
#  cbar=plt.colorbar(pad=0.01,ticks=np.linspace(0, 0.3, 4))
#  cbar.ax.tick_params(labelsize=font2['size']) 
#  cbar.set_label('$x|_{t=0}$',fontdict=font2)        
  plt.xlim(-0.5,4.5)
  plt.ylim(-7,7)
#  plt.xlabel('$x\ [\mu m]$',fontdict=font)
#  plt.ylabel('$y\ [\mu m]$',fontdict=font)
  plt.xticks(fontsize=0.001); 
  plt.yticks(fontsize=0.001);

  plt.subplots_adjust(left=0.15, bottom=0.15, right=0.92, top=0.98, wspace=0.01, hspace=0.01)
#      plt.title('Au51 at '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
  fig = plt.gcf()
  fig.set_size_inches(10, 6)
  fig.savefig('wrap_px_x_initial.png',format='png',dpi=160)
  plt.close("all")
  print('wrap_px_x_initial.png')
