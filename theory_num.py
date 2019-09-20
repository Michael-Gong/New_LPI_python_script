# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:49:31 2019

@author: kathe
"""

#!/public/home/users/bio001/tools/python-2.7.11/bin/python
#import sdf
import matplotlib
#matplotlib.use('agg')
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
  a0_th  = np.linspace(100,2000,200)
  a0_num = np.array([100,200,400,800,1200,1600,2000])

  tqe_num =np.array([650,280,100,20, 5,   2,   0])
  gqe_num =np.array([4.8,8.4,14, 0, 33,  0,   0])

  trr_num =np.array([1400,1000,600,100,50,10, 5 ])
  grr_num =np.array([3.2,4.1, 8.1, 0, 13,0,0])

  #I0 = a0**2*1.37e18
  #T0 = np.array([3,4,5,6,7,8])*1.0

  gg_th  = 4.4e3*a0_th**0.25/(5e-3*200.)**0.75
  tt_th  = 2.7e7*a0_th**(-3.)*a0_th/2.
  #t0 = 3.3/((alpha/(2*pi*rb))+(pi**2/(16*(rb**2))))
  plt.subplot(1,2,1)
  plt.scatter(a0_num, tqe_num,c='royalblue',marker='o',s=200, label='QED', edgecolors='black', linewidth='2',alpha=1,zorder=2)
  plt.scatter(a0_num, trr_num,c='royalblue',marker='^',s=200, label='Classic', edgecolors='black', linewidth='2',alpha=1,zorder=2)
  plt.plot(a0_th,tt_th*3.33,'--',c='k',linewidth=4, label='theory '+r'$\Delta t$',zorder=0)

  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('$a_0$',fontdict=font)
  plt.ylabel('t [fs]',fontdict=font)
  plt.xticks([0,500,1000,1500,2000],fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  #plt.yscale('log')
  plt.ylim(0,1500)
#  plt.xlim(1,7)
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.legend(loc='best',fontsize=20,framealpha=0.8)
#        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)


  plt.subplot(1,2,2)
  plt.scatter(a0_num, gqe_num,c='tomato',marker='o',s=200, label='QED', edgecolors='black', linewidth='2',alpha=1,zorder=2)
  plt.scatter(a0_num, grr_num,c='tomato',marker='^',s=200, label='Classic', edgecolors='black', linewidth='2',alpha=1,zorder=2)
  plt.plot(a0_th,gg_th*0.51e-3,'--',c='k',linewidth=4, label='theory '+r'$\gamma_m$',zorder=0)

  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('$a_0$',fontdict=font)
  plt.ylabel(r'$\varepsilon_e$'+' [GeV]',fontdict=font)
  plt.xticks([0,500,1000,1500,2000],fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  #plt.yscale('log')
#  plt.ylim(0,2750)
#  plt.xlim(1,7)
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.legend(loc='best',fontsize=20,framealpha=0.8)
#        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)


  plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.95,
            wspace=None, hspace=None)
  fig = plt.gcf()
  fig.set_size_inches(16, 7.6)
  #plt.show()
  fig.savefig('./theory_num.png',format='png',dpi=160)
  plt.close("all")
