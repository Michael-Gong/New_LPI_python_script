import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal

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
wavelength=     0.8e-6
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
 
font2 = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 18,  
        }  
 
font_size=25
font_size2=18

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

def function_plot_dN_dE(data,name,mass,linecolor):
#    print(data.keys())
    px = data['Particles/Px/'+name].data/(mass*m0*v0)
    py = data['Particles/Py/'+name].data/(mass*m0*v0)
    grid_y = data['Grid/Particles/'+name].data[1]/1e-6
    gg = (px**2+py**2+1.0)**0.5
    ek = (gg-1.)*mass*m0*v0**2/1.6e-13
    ww = data['Particles/Weight/'+name].data*4e-6
    theta = np.arctan2(py,px)*180.0/np.pi
    ek=ek[(px>0)&(abs(grid_y)<10)]
    ww=ww[(px>0)&(abs(grid_y)<10)]
#        ww_x_d = ww[ (grid_x>5) & (abs(grid_y)<3.) & (work_x>work_y)]
    dist_x1, den1 = pxpy_to_energy(ek,ww)
    plt.plot(dist_x1,den1,color=linecolor,linewidth=3,label=name)

def function_bar_work_frac(data,name,mass):
    work_x = data['Particles/Time_Integrated_Work_x/'+name].data  # value for gamma vector
    work_y = data['Particles/Time_Integrated_Work_y/'+name].data
    ek     = (work_x+work_y)*mass*0.51*1
    value_axisx = np.linspace(5,500,50)
    value_axisy = np.linspace(5,500,50)
    value_grid = np.linspace(0,500,51)
    value_total_x = np.zeros_like(value_axisy)
    value_total_y = np.zeros_like(value_axisy)
    value_num   = np.zeros_like(value_axisy)
    for i in range(np.size(value_axisy)):
        value_total_x[i] = np.sum(work_x[(value_grid[i]<=ek) & (value_grid[i+1]>ek)],0)
        value_total_y[i] = np.sum(work_y[(value_grid[i]<=ek) & (value_grid[i+1]>ek)],0)
        value_num[i] = np.size(work_y[(value_grid[i]<=ek) & (value_grid[i+1]>ek)])
        print('x-:',value_total_x[i]/(value_total_x[i]+value_total_y[i]),'; y-:',value_total_y[i]/(value_total_x[i]+value_total_y[i]))
    y_x = value_total_x/(value_total_x+value_total_y)
    #y_x[y_x > 1] = 1
    #y_y = 1-y_x
    y_y = value_total_y/(value_total_x+value_total_y)
    width=10
    pl=plt.bar(value_axisx, y_x*100, width, color='salmon',edgecolor='black',linewidth=1)
    pt=plt.bar(value_axisx, y_y*100, width, bottom=y_x*100, color='deepskyblue',edgecolor='black',linewidth=1)



if __name__ == "__main__":
  ######### Script code drawing figure ################
  n = 18
  mass=1836.
  name='proton' 
  plt.subplot(1,1,1)

  data = sdf.read('./PW_w010/'+str(n).zfill(4)+".sdf",dict=True)
  time=data['Header']['time']
  px = data['Particles/Px/'+name].data/(mass*m0*v0)
  py = data['Particles/Py/'+name].data/(mass*m0*v0)
  grid_y = data['Grid/Particles/'+name].data[1]/1e-6
  gg = (px**2+py**2+1.0)**0.5
  ek = (gg-1.)*mass*m0*v0**2/1.6e-13
  ww = data['Particles/Weight/'+name].data*4e-6
  theta = np.arctan2(py,px)*180.0/np.pi
  ek=ek[(px>0)&(abs(grid_y)<12)&(abs(theta)<30)]
  ww=ww[(px>0)&(abs(grid_y)<12)&(abs(theta)<30)]
  dist_x1, den1 = pxpy_to_energy(ek,ww)
  plt.plot(dist_x1,den1,color='crimson',linewidth=3,label='$\sigma_0=1\mu m$')

  data = sdf.read('./PW_w020/'+str(n).zfill(4)+".sdf",dict=True)
  time=data['Header']['time']
  px = data['Particles/Px/'+name].data/(mass*m0*v0)
  py = data['Particles/Py/'+name].data/(mass*m0*v0)
  grid_y = data['Grid/Particles/'+name].data[1]/1e-6
  gg = (px**2+py**2+1.0)**0.5
  ek = (gg-1.)*mass*m0*v0**2/1.6e-13
  ww = data['Particles/Weight/'+name].data*4e-6
  theta = np.arctan2(py,px)*180.0/np.pi
  ek=ek[(px>0)&(abs(grid_y)<12)&(abs(theta)<30)]
  ww=ww[(px>0)&(abs(grid_y)<12)&(abs(theta)<30)]
  dist_x1, den1 = pxpy_to_energy(ek,ww)
  plt.plot(dist_x1,den1,color='seagreen',linewidth=3,label='$\sigma_0=2\mu m$')

  data = sdf.read('./PW_w040/'+str(n).zfill(4)+".sdf",dict=True)
  time=data['Header']['time']
  px = data['Particles/Px/'+name].data/(mass*m0*v0)
  py = data['Particles/Py/'+name].data/(mass*m0*v0)
  grid_y = data['Grid/Particles/'+name].data[1]/1e-6
  gg = (px**2+py**2+1.0)**0.5
  ek = (gg-1.)*mass*m0*v0**2/1.6e-13
  ww = data['Particles/Weight/'+name].data*4e-6
  theta = np.arctan2(py,px)*180.0/np.pi
  ek=ek[(px>0)&(abs(grid_y)<12)&(abs(theta)<30)]
  ww=ww[(px>0)&(abs(grid_y)<12)&(abs(theta)<30)]
  dist_x1, den1 = pxpy_to_energy(ek,ww)
  plt.plot(dist_x1,den1,color='mediumblue',linewidth=3,label='$\sigma_0=4\mu m$')

  data = sdf.read('./PW_w080/'+str(n).zfill(4)+".sdf",dict=True)
  time=data['Header']['time']
  px = data['Particles/Px/'+name].data/(mass*m0*v0)
  py = data['Particles/Py/'+name].data/(mass*m0*v0)
  grid_y = data['Grid/Particles/'+name].data[1]/1e-6
  gg = (px**2+py**2+1.0)**0.5
  ek = (gg-1.)*mass*m0*v0**2/1.6e-13
  ww = data['Particles/Weight/'+name].data*4e-6
  theta = np.arctan2(py,px)*180.0/np.pi
  ek=ek[(px>0)&(abs(grid_y)<12)&(abs(theta)<30)]
  ww=ww[(px>0)&(abs(grid_y)<12)&(abs(theta)<30)]
  dist_x1, den1 = pxpy_to_energy(ek,ww)
  plt.plot(dist_x1,den1,color='black',linewidth=3,label='$\sigma_0=8\mu m$')

  plt.xlabel(r'$\varepsilon_p$'+' [MeV]',fontdict=font)
  plt.ylabel('dN/dE [A.U.]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.yscale('log')
  plt.xlim(0,390)
#    plt.ylim(1e1,1e5)
  plt.legend(loc='best',fontsize=font_size2,framealpha=1)
  plt.title('t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
  plt.grid(which='major',linestyle=':', linewidth=1, color='k')

  plt.subplots_adjust(left=0.15, bottom=0.18, right=0.98, top=0.95, wspace=0.011, hspace=0.16)
  fig = plt.gcf()
  fig.set_size_inches(10, 5.)
  fig.savefig('wrap_dN_dE.png',format='png',dpi=160)
  plt.close("all")
  print('wrap_dN_dE.png')

