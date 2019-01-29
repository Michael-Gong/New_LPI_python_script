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
font_size = 25
######### Parameter you should set ###########

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
    lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


if __name__ == "__main__":
  from_path = './two/'
  to_path   = './two/'
#  nsteps      = 79580 #sum(1 for line in open(from_path+'x_0000.txt'))/part_number

  xx  =  np.load(to_path+'xx2d.npy')
  yy  =  np.load(to_path+'yy2d.npy')

  ex  =  np.zeros_like(xx)
  ey  =  np.zeros_like(xx)
  bz  =  np.zeros_like(xx)

  part_number = len(xx[:,0])
  nsteps      = len(xx[0,:])
  for n in range(part_number):
      for i in range(nsteps):
          data_e  = sdf.read(from_path+'e_fields'+str(i).zfill(4)+'.sdf',dict=True)
          data_b  = sdf.read(from_path+'b_fields'+str(i).zfill(4)+'.sdf',dict=True)
          grid_x0 = data_b['Grid/Grid'].data[0]/1.0e-6
          grid_y0 = data_b['Grid/Grid'].data[1]/1.0e-6
          grid_x  = 0.5*(grid_x0[:-1]+grid_x0[1:])   
          grid_y  = 0.5*(grid_y0[:-1]+grid_y0[1:])   
 
          tx  = xx[n,i]
          ty  = yy[n,i]
          ii = np.max(np.where( xx[n,i] > grid_x))
          jj = np.max(np.where( yy[n,i] > grid_y))
          d_xy = (grid_x[-1]-grid_x[-2])*(grid_y[-1]-grid_y[-2])
         
          data  = data_e['Electric Field/Ey'].data/exunit
          ey[n,i] = data[ ii,  jj ]*(grid_x[ii+1]-tx)*(grid_y[jj+1]-ty)/d_xy\
                   +data[ii+1, jj ]*(tx - grid_x[ii])*(grid_y[jj+1]-ty)/d_xy\
                   +data[ ii, jj+1]*(grid_x[ii+1]-tx)*(ty - grid_y[jj])/d_xy\
                   +data[ii+1,jj+1]*(tx - grid_x[ii])*(ty - grid_y[jj])/d_xy
    
          data  = data_e['Electric Field/Ex'].data/exunit
          ex[n,i] = data[ ii,  jj ]*(grid_x[ii+1]-tx)*(grid_y[jj+1]-ty)/d_xy\
                   +data[ii+1, jj ]*(tx - grid_x[ii])*(grid_y[jj+1]-ty)/d_xy\
                   +data[ ii, jj+1]*(grid_x[ii+1]-tx)*(ty - grid_y[jj])/d_xy\
                   +data[ii+1,jj+1]*(tx - grid_x[ii])*(ty - grid_y[jj])/d_xy
    
          data  = data_b['Magnetic Field/Bz'].data/bxunit
          bz[n,i] = data[ ii,  jj ]*(grid_x[ii+1]-tx)*(grid_y[jj+1]-ty)/d_xy\
                   +data[ii+1, jj ]*(tx - grid_x[ii])*(grid_y[jj+1]-ty)/d_xy\
                   +data[ ii, jj+1]*(grid_x[ii+1]-tx)*(ty - grid_y[jj])/d_xy\
                   +data[ii+1,jj+1]*(tx - grid_x[ii])*(ty - grid_y[jj])/d_xy
     
          print('finished! '+str(n)+' :'+str(i).zfill(4))


  np.save(to_path+'ey2d',ey)
  np.save(to_path+'ex2d',ex)
  np.save(to_path+'bz2d',bz)
