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
        'size'   : 20,  
       }  

font2 = {'family' : 'monospace',  
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 20,  
       }

start   =  0
stop    =  850
step    =  1

from_path = './two/'
to_path  = from_path 

ey_vph    = np.zeros([4000,stop-start+1])
grid_time = np.zeros(stop-start+1)
grid_x    = np.zeros(4000)

for n in range(start,stop+step,step):
    data = sdf.read(from_path+'e_fields'+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid'].data[0]/1.0e-6
    name = 'ey'
    eexx = data['Electric Field/'+str.capitalize(name)].data/exunit
    #n3d=len(eexx[0,0,:])
    #ex = np.sum(eexx[:,:,(n3d//2-1):(n3d//2+1)],2)/2
    #print(np.shape(ex))
    ex = eexx
    n3d=len(ex[0,:])
    ex = np.sum(ex[:,(n3d//2-1):(n3d//2+1)],1)/2
    print(np.shape(ex))
    ey_vph[:,n]=ex
    grid_time[n] = time
    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')

grid_x = 0.5*(x[:-1]+x[1:])
np.save(to_path+'ey_vph',ey_vph)
np.save(to_path+'grid_time',grid_time)
np.save(to_path+'grid_x',grid_x)


