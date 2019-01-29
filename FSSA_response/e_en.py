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
        'size'   : 20,  
        }  


def processplot(n): 
  ######### Parameter you should set ###########
  #start   =  23  # start time
  #stop    =  30  # end time
  #step    =  1  # the interval or step
  
  #youwant = ['E_u_1_density','E_u_density','Ion_density','Ion_1_density']
  #youwant = ['Electron_density','E_1_density']
  #youwant =  ['Electron_ekbar']#['Electron_density','E_1_density','ey','ex','ey_averaged','bz','bz_averaged','ex_averaged','Electron_en']
  youwant =  ['Electron_en']
  #youwant =  ['ey','ex','ey_averaged','bz','bz_averaged','Electron_ekbar','Electron_density','Ion_density','Ion_ekbar','Carbon_density','Carbon_ekbar','ex_averaged']
  #youwant.append('Ion_ekbar')
  #youwant.append('positron_ekbar')
  #youwant.append('electron_en')
  #youwant.append('photon_en')
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px...
  
  from_path = './one_0.25pi/'
  to_path=from_path
  
  ######### Script code drawing figure ################
  #for n in range(start,stop+step,step):
  #### header data ####
  #data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
  #header=data['Header']
  #time=header['time']
  #x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  #y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  #X, Y = np.meshgrid(x, y)
  
  data = sdf.read(from_path+'dist'+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time=header['time']
  name = 'Electron_en'
  den = data['dist_fn/en/'+name[0:-3]].data*6.4e-6
  dist_x  = data['Grid/en/'+name[0:-3]].data[0]/(q0*1.0e6)
  plt.plot(dist_x,den/(dist_x[1]-dist_x[0]),'-',color='salmon',linewidth=3)
  plt.text(100,1e13,'>1MeV  [nC]'+str(np.sum(den[dist_x>1])*1.6e-19/1e-9))
  plt.text(100,3.2e12,'>10MeV [nC]'+str(np.sum(den[dist_x>10])*1.6e-19/1e-9))
  plt.text(100,1e12,'>100MeV[nC]'+str(np.sum(den[dist_x>100])*1.6e-19/1e-9))
  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('Energy [MeV]',fontdict=font)
  plt.ylabel('dN/dE [MeV$^-1$]',fontdict=font)
  plt.xticks(fontsize=20); plt.yticks(fontsize=20);
  plt.yscale('log')
  plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  fig = plt.gcf()
  fig.set_size_inches(12, 7)
  fig.savefig(to_path+'e_en_'+str(n).zfill(4)+'.png',format='png',dpi=100)
  plt.close("all")
  print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
  return 0

if __name__ == '__main__':
  start   =  0 # start time
  stop    =  32  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=1)
  results = pool.map(processplot,inputs)
  print(results)
