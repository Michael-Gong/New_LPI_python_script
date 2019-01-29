
#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
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
  wavelength=     1.0e-6
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
          'size'   : 20,  
          }

  from_path_0=['./Data_l0600/' ,'./Data_l0900/','./Data_l1200/','./Data_l1380/','./Data_l1500/','./Data_l1640/','./Data_l1800/','./Data_l2000/','./Data_l2240/','./Data_l2560/']


for i in range(np.size(from_path_0)):
  from_path = from_path_0[i]
  for n in range(14):
    data_1 = sdf.read(from_path+'b_fields'+str(n).zfill(4)+'.sdf',dict=True)
    data_2 = sdf.read(from_path+'e_fields'+str(n).zfill(4)+'.sdf',dict=True)
   
    bz = data_1['Magnetic Field/Bz'].data/bxunit
    ey = data_2['Electric Field/Ey'].data/exunit

    x,y = data_2['Grid/Grid_mid'].data
    X,Y = np.meshgrid(x,y,indexing='ij')
 

    forward  = np.sum((ey+bz)**2/4.0)
    backward = np.sum((ey-bz)**2/4.0)

    if n==0:
        print('Ratio for '+from_path)
        energy0 = forward
     
    print('For '+str(n).zfill(4)+': forward '+str(forward)+'='+str(forward/energy0*100)+'%'+'; backward '+str(backward)+'='+str(backward/energy0*100)+'%')

