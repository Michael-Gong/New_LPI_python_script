#!/public/home/users/bio001/tools/python-2.7.11/bin/python
#import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
  
if __name__ == "__main__":
  print ('This is main of module "test2d.py"')
  ######## Constant defined here ########
  pi        =     3.1415926535897932384626
  q0        =     1.602176565e-19 # C
  m0        =     9.10938291e-31  # kg
  v0        =     2.99792458e8    # m/s^2
  kb        =     1.3806488e-23   # J/K
  mu0       =     4.0e-7*np.pi       # N/A^2
  epsilon0  =     8.8541878176203899e-12 # F/m
  h_planck  =     6.62606957e-34  # J s
  wavelength=     1.06e-6
  frequency =     v0*2*pi/wavelength
  
  exunit    =     m0*v0*frequency/q0
  bxunit    =     m0*frequency/q0
  denunit    =     frequency**2*epsilon0*m0/q0**2
  jalf      =     4*np.pi*epsilon0*m0*v0**3/q0/wavelength**2
  print('electric field unit: '+str(exunit))
  print('magnetic field unit: '+str(bxunit))
  print('density unit nc: '+str(denunit))
  
  font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 25,  
          }  
  font_size=25

##below is for norm colorbar
  class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y)) 
##end for norm colorbar####


  def pxpy_to_energy(gamma, weight):
      binsize = 500
      en_grid = np.linspace(0.5,499.5,binsize)
      en_bin  = np.linspace(0,500.0,binsize+1)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = np.sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/(500.0/binsize)
      return (en_grid, en_value)

  def theta_formula(a0, tau0, c1, c2):
      w0 = 5e-6
      return (4*np.pi*v0)*m0*tau0*a0/(mu0*q0**2*wavelength*w0**2*20*denunit*(1.+1836.*c2*a0**(-0.5)))*c1

  to_path='./'
  a0 = np.array([2.0, 2.83, 4.00, 5.66, 8.00])
  I0 = a0**2*1.37e18
  T0 = np.array([3,4,5,6,7,8])*1.0

  T0_0 = np.array([-5.79e-1, -8.03e-1, -1.35, -2.34, -3.66])
  T0_1 = np.array([-6.74e-1, -9.65e-1, -1.63, -2.87, -4.23])
  T0_2 = np.array([-8.16e-1, -1.23,    -2.0 , -3.33, -4.73])
  T0_3 = np.array([-9.84e-1, -1.51,    -2.46, -3.74, -5.44])
  T0_4 = np.array([-1.19,    -1.78,    -2.73, -4.36, -5.95])
  T0_5 = np.array([-1.38,    -2.11,    -2.98, -4.78, -6.88])

#  plt.scatter(a0,-T0_0,c='silver',marker='o',s=100, label='$t=3T_0$', edgecolors='black', linewidth='2',alpha=1,zorder=2)
  plt.scatter(a0,-T0_1,c='blueviolet',marker='o',s=200, label='$t=4T_0$', edgecolors='black', linewidth='2',alpha=1,zorder=2)
#  plt.scatter(a0,-T0_2,c='deepskyblue',marker='o',s=100, label='$t=5T_0$', edgecolors='black', linewidth='2',alpha=1,zorder=2)
  plt.scatter(a0,-T0_3,c='limegreen',marker='o',s=200, label='$t=6T_0$', edgecolors='black', linewidth='2',alpha=1,zorder=2)
#  plt.scatter(a0,-T0_4,c='orange',marker='o',s=100, label='$t=7T_0$', edgecolors='black', linewidth='2',alpha=1,zorder=2)
  plt.scatter(a0,-T0_5,c='tomato',marker='o',s=200, label='$t=8T_0$', edgecolors='black', linewidth='2',alpha=1,zorder=2)

  a0=np.linspace(1,10,101)
  plt.plot(a0, theta_formula(a0, 6*3.3e-15, 0.42, 0.002)/1e-3, '-k',linewidth=4, label='$C_1=0.42,C_2=0.02$',zorder=0)
  plt.plot(a0, theta_formula(a0, 4*3.3e-15, 0.42, 0.002)/1e-3, '--k',linewidth=4, label='$C_1=0.42,C_2=0.02$',zorder=0)
  plt.plot(a0, theta_formula(a0, 8*3.3e-15, 0.42, 0.002)/1e-3, ':k',linewidth=4, label='$C_1=0.42,C_2=0.02$',zorder=0)
#  plt.plot(a0, theta_formula(a0, 3*3.3e-15, 0.35, 0.0015)/1e-3, '--k',linewidth=4, label='$C_1=0.19,C_2=0.1$',zorder=0)
#  plt.plot(a0, theta_formula(a0, 5*3.3e-15, 0.35, 0.0015)/1e-3, '--k',linewidth=4, label='$C_1=0.19,C_2=0.1$',zorder=0)

  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('$a_0$',fontdict=font)
  plt.ylabel(r'$\Delta\Theta$'+' [mmrad]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  #plt.yscale('log')
  #plt.ylim(2e7,8e9)
  plt.xlim(1,9)
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)

  plt.legend(loc='best',fontsize=20,framealpha=0.8)
  plt.subplots_adjust(left=None, bottom=0.15, right=0.95, top=0.95,
            wspace=None, hspace=None)
#        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
  fig = plt.gcf()
  fig.set_size_inches(8, 6.8)
  fig.savefig('./crater_scan.png',format='png',dpi=160)
  plt.close("all")
