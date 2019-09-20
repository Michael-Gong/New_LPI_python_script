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

upper = matplotlib.cm.YlOrBr(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_or = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

upper = matplotlib.cm.autumn(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_autumn = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


def pxpy_to_energy(gamma, weight):
    binsize = 400
    en_grid = np.linspace(20,15980,400)
    en_bin  = np.linspace(0,16000.0,401)
    en_value = np.zeros_like(en_grid)
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
    return (en_grid, en_value)


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
  gg_no_0 = (px0**2+py0**2+1)**0.5

  gg_no = convert_g_m(gg_no_0)



  from_path = './part_01_rr/'
  t0  = np.loadtxt(from_path+'t_0000.txt')/2/np.pi
  px0 = np.loadtxt(from_path+'px_0000.txt')
  py0 = np.loadtxt(from_path+'py_0000.txt')
  t0  = np.reshape(t0,(part_number,nsteps))
  px0 = np.reshape(px0,(part_number,nsteps))
  py0 = np.reshape(py0,(part_number,nsteps))
  gg_rr_0 = (px0**2+py0**2+1)**0.5

  gg_rr = convert_g_m(gg_rr_0)


  part_number=5000
  from_path = './part_01_qe/'
  nsteps      = int(sum(1 for line in open(from_path+'t_0000.txt'))/part_number)
  t1  = np.loadtxt(from_path+'t_0000.txt')/2/np.pi
  px0 = np.loadtxt(from_path+'px_0000.txt')
  py0 = np.loadtxt(from_path+'py_0000.txt')
  t1  = np.reshape(t1,(part_number,nsteps))
  px0 = np.reshape(px0,(part_number,nsteps))
  py0 = np.reshape(py0,(part_number,nsteps))
  gg_qe_0 = (px0**2+py0**2+1)**0.5
  ww  = np.zeros_like(gg_qe_0)+1

  gg_qe = convert_g_m(gg_qe_0)


  axis_time = t1[0,:]
  print(np.shape(axis_time))
  axis_en   = np.linspace(20,15980,400)
  data_qe   = np.zeros([np.size(axis_time),np.size(axis_en)])

  for i in range(nsteps): 
      axis_en, data_qe[i,:] = pxpy_to_energy(gg_qe[:,i],ww[:,i]) 

  plt.subplot(1,1,1)
#  for n in range(25,75,10):
  x,y=np.meshgrid(axis_time,axis_en)
  levels = np.linspace(0, 100, 41)
  plt.pcolormesh(x, y/1e3, data_qe.T/5000.0*100.0, norm=colors.Normalize(vmin=0, vmax=5), cmap=mycolor_autumn)
#plt.pcolormesh(x, y, ex.T, norm=mpl.colors.Normalize(vmin=0,vmax=100,clip=True), cmap=cm.cubehelix_r)
#  plt.axis([x.min(), x.max(), y.min(), y.max()])
#### manifesting colorbar, changing label and axis properties ####
  cbar=plt.colorbar(pad=0.01)#ticks=[np.min(ex), -eee/2, 0, eee/2, np.min()])
  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
  cbar.set_label('h [%]',fontdict=font)

  n=0
  #plt.scatter((t0-x0)[n,:], (gg0-px0)[n,:], c=np.zeros_like(px0[n,:])+py0[n,0], norm=colors.Normalize(vmin=50,vmax=150), s=5, cmap='rainbow', edgecolors='None', alpha=1)
  plt.plot(t0[n,:], gg_no[n,:]/1000,'-',color='k',linewidth=3)
  plt.plot(t0[n,:], gg_rr[n,:]/1000,'-',color='b',linewidth=3)
#  for n in range(10000):
#      plt.plot(t1[n,:], gg_qe[n,:]/1000,'-',color='r',linewidth=0.01, zorder=0)
 # plt.plot(t0[n,:], gg_no_0[n,:],'-',color='k',linewidth=3, zorder=0)
 # plt.plot(t0[n,:], gg_rr_0[n,:],'-',color='b',linewidth=3, zorder=0)
 #   plt.legend(loc='upper right')

 #   plt.legend(loc='upper right')
  plt.xlim(0,320)
  plt.ylim(0,15)
  plt.xlabel('t $[T_0]$',fontdict=font)
  plt.ylabel(r'$\varepsilon_e\ [10^3m_ec^2]$',fontdict=font)
#  plt.xticks([100,200,300],fontsize=font_size)
#  plt.yticks([5,10,15],fontsize=font_size)
  plt.xticks(fontsize=font_size)
  plt.yticks(fontsize=font_size)
#  plt.title('t='+str(round(t0[0,i],0))+' $T_0$',fontdict=font)
  #plt.text(-100,650,' t = 400 fs',fontdict=font)

  plt.subplots_adjust(left=0.15, bottom=0.16, right=0.98, top=0.98,
                wspace=None, hspace=None)

  plt.show()
  #lt.figure(figsize=(100,100))
  fig = plt.gcf()
  fig.set_size_inches(10, 8.5)
  fig.savefig('max_5000_con.png',format='png',dpi=160)
  plt.close("all")
