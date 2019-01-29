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

font_size = 25

if __name__ == '__main__':
  ######### Parameter you should set ###########
  #start   =  23  # start time
  #stop    =  30  # end time
  #step    =  1  # the interval or step
  
  #youwant = ['E_u_1_density','E_u_density','Ion_density','Ion_1_density']
  #youwant = ['Electron_density','E_1_density']
  #youwant =  ['Electron_ekbar']#['Electron_density','E_1_density','ey','ex','ey_averaged','bz','bz_averaged','ex_averaged','Electron_en']


  from_path = './one_double/'

  plt.subplot2grid((1, 2), (0, 0))
  n=13
  data = sdf.read(from_path+'current'+str(n).zfill(4)+".sdf",dict=True)
  x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  print('ok')
  y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  x1 = np.linspace(-5,75,600)
  J1 = 0.3*(2*np.pi)**2*3.2**2 
  J2 = 0.056*(2*np.pi)**2*3.2**2 
  y1 = np.zeros_like(x1)
  y1[(x1>10)&(x1<45)] = J1
  y2 = np.zeros_like(x1)+J2
    
  plt.plot(x1, y2, '--r', linewidth=3,label=r'$|J_x|=\alpha^*(\frac{2\pi r}{\lambda})^2J_A,\ \alpha^*=0.056$',zorder=0)
#    ax.plot(x1, y1, '--b', linewidth=3,label=r'$-J_x=\alpha(2\pi r/\lambda)^2J_A,\ \ \ \alpha=0.30$',zorder=1)
  data = sdf.read(from_path+'current'+str(n).zfill(4)+".sdf",dict=True)
  ex = data['Current/Jx_averaged'].data
  Jx = np.sum(ex[:,abs(y)<1.6],1)*(1e-6/30)*6.4e-6/17e3
  plt.plot(x, -Jx, '-b', linewidth=3,label='2D-PIC simulation')
  plt.xlim(-5,55)
  plt.xlabel('X [$\mu m$]',fontdict=font)
  plt.ylabel('$-J_x$ [$J_A$]',fontdict=font)
  plt.xticks(fontsize=font_size) 
  plt.yticks(fontsize=font_size) 
  plt.legend(loc='upper left',fontsize=20,framealpha=0.0)
  plt.ylim(0,120)

  plt.subplot2grid((1, 2), (0, 1))
  
  youwant =  ['Electron_en']
  #youwant =  ['ey','ex','ey_averaged','bz','bz_averaged','Electron_ekbar','Electron_density','Ion_density','Ion_ekbar','Carbon_density','Carbon_ekbar','ex_averaged']
  #youwant.append('Ion_ekbar')
  #youwant.append('positron_ekbar')
  #youwant.append('electron_en')
  #youwant.append('photon_en')
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px...
  
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
  n = 8  
  data = sdf.read(from_path+'dist'+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time1=header['time']
  name = 'Electron_en'
  den1 = data['dist_fn/en/'+name[0:-3]].data*6.4e-6/0.95
  dist_x1  = data['Grid/en/'+name[0:-3]].data[0]/(q0*1.0e6)
  n = 13 
  data = sdf.read(from_path+'dist'+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time2=header['time']
  name = 'Electron_en'
  den2 = data['dist_fn/en/'+name[0:-3]].data*6.4e-6/0.95
  dist_x2  = data['Grid/en/'+name[0:-3]].data[0]/(q0*1.0e6)
  n = 16  
  data = sdf.read(from_path+'dist'+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time3=header['time']
  name = 'Electron_en'
  den3 = data['dist_fn/en/'+name[0:-3]].data*6.4e-6/0.95
  dist_x3  = data['Grid/en/'+name[0:-3]].data[0]/(q0*1.0e6)

  plt.plot(dist_x1,den1/(dist_x1[1]-dist_x1[0]),'-',color='royalblue',linewidth=3,label='t = '+str(round(time1/1.0e-15,0))+' fs')
  plt.plot(dist_x2,den2/(dist_x2[1]-dist_x2[0]),'-',color='limegreen',linewidth=3,label='t = '+str(round(time2/1.0e-15,0))+' fs')
  plt.plot(dist_x3,den3/(dist_x3[1]-dist_x3[0]),'-',color='tomato',linewidth=3,label='t = '+str(round(time3/1.0e-15,0))+' fs')

#  plt.text(100,1e13,'>1MeV  [nC]'+str(np.sum(den[dist_x>1])*1.6e-19/1e-9))
#  plt.text(100,3.2e12,'>10MeV [nC]'+str(np.sum(den[dist_x>10])*1.6e-19/1e-9))
#  plt.text(100,1e12,'>100MeV[nC]'+str(np.sum(den[dist_x>100])*1.6e-19/1e-9))
  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('Energy [MeV]',fontdict=font)
  plt.ylabel('dN/dE [MeV$^-1$]',fontdict=font)
  plt.xticks([0,200,400,600,800],fontsize=font_size); 
  plt.yticks(np.logspace(6,12,4),fontsize=font_size);
  plt.yscale('log')
#  plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
#  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.legend(loc='best',fontsize=20,framealpha=0.0)
  plt.xlim(0,850); plt.ylim(5e5,1e12)



  plt.subplots_adjust(left=0.08, bottom=0.14, right=0.99, top=0.95, wspace=0.25, hspace=0.025)


  fig = plt.gcf()
  fig.set_size_inches(15, 7)
  fig.savefig(to_path+'wrap_en.png',format='png',dpi=160)
  plt.close("all")
  print('finised %')

