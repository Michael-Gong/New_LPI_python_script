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


def convert_g_m(gg):
    gg_m = np.zeros_like(gg)
    for i in range(np.size(gg[:,0])):  
        for j in range(np.size(gg[0,:])):
            gg_m[i,j] = np.max(gg[i,:(j+1)])
    return gg_m

if __name__ == "__main__":
  part_number = 1
  from_path = './part_01_no/'
  nsteps      = int(sum(1 for line in open(from_path+'t_0000.txt'))/part_number)

  from_path = './part_01_no/'
  t0  = np.loadtxt(from_path+'t_0000.txt')/2/np.pi
  px0 = np.loadtxt(from_path+'px_0000.txt')
  py0 = np.loadtxt(from_path+'py_0000.txt')
  t0  = np.reshape(t0,(part_number,nsteps))
  px0 = np.reshape(px0,(part_number,nsteps))
  py0 = np.reshape(py0,(part_number,nsteps))
  gg0 = (px0**2+py0**2+1)**0.5

  gg_no = convert_g_m(gg0)



  from_path = './part_01_rr/'
  t0  = np.loadtxt(from_path+'t_0000.txt')/2/np.pi
  px0 = np.loadtxt(from_path+'px_0000.txt')
  py0 = np.loadtxt(from_path+'py_0000.txt')
  t0  = np.reshape(t0,(part_number,nsteps))
  px0 = np.reshape(px0,(part_number,nsteps))
  py0 = np.reshape(py0,(part_number,nsteps))
  gg0 = (px0**2+py0**2+1)**0.5

  gg_rr = convert_g_m(gg0)


  part_number=5000
  from_path = './part_01_qe/'
  px0 = np.loadtxt(from_path+'px_0000.txt')
  py0 = np.loadtxt(from_path+'py_0000.txt')
  px0 = np.reshape(px0,(part_number,nsteps))
  py0 = np.reshape(py0,(part_number,nsteps))
  gg0 = (px0**2+py0**2+1)**0.5

  gg_qe = convert_g_m(gg0)



  plt.subplot(1,1,1)
#  for n in range(25,75,10):
  n=0
  #plt.scatter((t0-x0)[n,:], (gg0-px0)[n,:], c=np.zeros_like(px0[n,:])+py0[n,0], norm=colors.Normalize(vmin=50,vmax=150), s=5, cmap='rainbow', edgecolors='None', alpha=1)
  plt.plot(t0[n,:], gg_no[n,:],':',color='k',linewidth=3, zorder=0)
  plt.plot(t0[n,:], gg_rr[n,:],':',color='b',linewidth=3, zorder=0)
 #   plt.legend(loc='upper right')

 #   plt.legend(loc='upper right')
 # plt.xlim(-0.2,30.2)
 # plt.ylim(0,105)
  plt.xlabel('t',fontdict=font)
  plt.ylabel(r'$\varepsilon_e [m_ec^2]$',fontdict=font)
  plt.xticks(fontsize=font_size)
  plt.yticks(fontsize=font_size);
#  plt.title('t='+str(round(t0[0,i],0))+' $T_0$',fontdict=font)
  #plt.text(-100,650,' t = 400 fs',fontdict=font)

  plt.subplots_adjust(left=0.15, bottom=0.16, right=0.98, top=0.98,
                wspace=None, hspace=None)

  plt.show()
  #lt.figure(figsize=(100,100))
  fig = plt.gcf()
  fig.set_size_inches(10, 8.5)
  fig.savefig('max_gg_t.png',format='png',dpi=160)
  plt.close("all")
