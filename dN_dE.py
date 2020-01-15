import sdf
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
 
  font_size=20 
  font_size2=14 
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
      binsize = 100
      en_grid = np.linspace(5,995,100)
      en_bin  = np.linspace(0,1000,101)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
      return (en_grid, en_value)

  from_path='./PW_w010/'
  to_path='./PW_w010_fig/'

  ######### Parameter you should set ###########
  start   =  1  # start time
  stop    =  19  # end time
  step    =  1  # the interval or step
  
  ######### Script code drawing figure ################
  for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    print('reading',from_path+str(n).zfill(4)+".sdf")
    
    px = data['Particles/Px/carbon'].data/(1836.*12.*m0*v0)
    py = data['Particles/Py/carbon'].data/(1836.*12.*m0*v0)
    grid_y = data['Grid/Particles/carbon'].data[1]/1e-6
    gg = (px**2+py**2+1.0)**0.5
    ek = (gg-1.)*1836.*12.*m0*v0**2/1.6e-13
    ww = data['Particles/Weight/carbon'].data*4e-6
    theta = np.arctan2(py,px)*180.0/np.pi
    ek=ek[(px>0)&(abs(grid_y)<10)]
    ww=ww[(px>0)&(abs(grid_y)<10)]
#        ww_x_d = ww[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
    dist_x1, den1 = pxpy_to_energy(ek,ww)
    plt.plot(dist_x1,den1,color='limegreen',linewidth=3,label='carbon')

    px = data['Particles/Px/electron'].data/(m0*v0)
    py = data['Particles/Py/electron'].data/(m0*v0)
    grid_y = data['Grid/Particles/electron'].data[1]/1e-6
    gg = (px**2+py**2+1.0)**0.5
    ek = (gg-1.)*m0*v0**2/1.6e-13
    ww = data['Particles/Weight/electron'].data*4e-6
    theta = np.arctan2(py,px)*180.0/np.pi
    ek=ek[(px>0)&(abs(grid_y)<10)]
    ww=ww[(px>0)&(abs(grid_y)<10)]
#        ww_x_d = ww[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
    dist_x1, den1 = pxpy_to_energy(ek,ww)
    plt.plot(dist_x1,den1,color='blue',linewidth=3,label='electron')

    px = data['Particles/Px/proton'].data/(1836.*m0*v0)
    py = data['Particles/Py/proton'].data/(1836.*m0*v0)
    grid_y = data['Grid/Particles/proton'].data[1]/1e-6
    gg = (px**2+py**2+1.0)**0.5
    ek = (gg-1.)*1836.*m0*v0**2/1.6e-13
    ww = data['Particles/Weight/proton'].data*4e-6
    theta = np.arctan2(py,px)*180.0/np.pi
    ek=ek[(px>0)&(abs(grid_y)<10)]
    ww=ww[(px>0)&(abs(grid_y)<10)]
#        ww_x_d = ww[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
    dist_x1, den1 = pxpy_to_energy(ek,ww)
    plt.plot(dist_x1,den1,color='red',linewidth=3,label='proton')

    plt.xlabel('Energy [MeV]',fontdict=font)
    plt.ylabel('dN/dE [A.U.]',fontdict=font)
    plt.xticks(fontsize=font_size); 
    plt.yticks(fontsize=font_size);
    plt.yscale('log')
    plt.xlim(0,500)
#    plt.ylim(1e1,1e5)
    plt.legend(loc='best',fontsize=4,framealpha=0.5)
    plt.title('t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    fig = plt.gcf()
    fig.set_size_inches(10, 7.)
    fig.savefig(to_path+'en_'+str(n).zfill(4)+'.png',format='png',dpi=320)
    plt.close("all")
    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
