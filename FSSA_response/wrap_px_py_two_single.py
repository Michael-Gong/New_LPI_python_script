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
  to_path = './two/'
#  nsteps      = 79580 #sum(1 for line in open(from_path+'x_0000.txt'))/part_number

  px  =  np.load(to_path+'px2d.npy')
  py  =  np.load(to_path+'py2d.npy')
  xx  =  np.load(to_path+'xx2d.npy')
  yy  =  np.load(to_path+'yy2d.npy')
  work_x = np.load(to_path+'workx2d.npy')
  work_y = np.load(to_path+'worky2d.npy')
  #gg  =  (px**2+py**2+pz**2+1.0)**0.5
  gg  =  (px**2+py**2+1.0)**0.5

  part_number = len(px[:,0])
  nsteps      = len(px[0,:])
  for i in range(np.size(xx[:,0])):
      n=7
      plt.subplot(2,1,1)
#      plt.plot(px[n,:], py[n,:], color='k', linewidth=3, linestyle='-')
      plt.plot(px[i,:], py[i,:], color='k', linewidth=3, linestyle='-')
      plt.scatter(px[i,:], py[i,:], c=gg[n,:], norm=colors.Normalize(vmin=0,vmax=1200), s=10, cmap='jet', edgecolors='None', alpha=1,zorder=3)
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.xlabel('$p_x$ [$m_ec$]',fontdict=font)
      plt.ylabel('$p_y$ [$m_ec$]',fontdict=font)
      plt.xticks([0,400,800,1200],fontsize=font_size); plt.yticks([-200,0,200],fontsize=font_size);
      #plt.xscale('log')
      plt.xlim(-50,1400)
      plt.ylim(-380,380)
      plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)

      plt.subplot(2,1,2)
      line_y1 = np.linspace(-12,-3.2,1001)
      line_x1 = np.zeros_like(line_y1)
    
      line_y2 = np.linspace(3.2,12,1001)
      line_x2 = np.zeros_like(line_y2)
    
      line_x3 = np.linspace(0,100.0,1001)
      line_y3 = np.zeros_like(line_x3)+3.2
    
      line_x4 = np.linspace(0,100.0,1001)
      line_y4 = np.zeros_like(line_x4)-3.2
      plt.plot(xx[i,:], yy[i,:], color='k', linewidth=3, linestyle='-')
      plt.scatter(xx[i,:], yy[i,:], c=gg[i,:], norm=colors.Normalize(vmin=0,vmax=1200), s=10, cmap='jet', edgecolors='None', alpha=1,zorder=3)
#  cbar=plt.colorbar( ticks=np.linspace(0, 500, 3) ,pad=0.005)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('$\gamma$',fontdict=font)
      plt.plot(line_x1,line_y1,linewidth=3,linestyle=':',color='k')
      plt.plot(line_x2,line_y2,linewidth=3,linestyle=':',color='k')
      plt.plot(line_x3,line_y3,linewidth=3,linestyle=':',color='k')
      plt.plot(line_x4,line_y4,linewidth=3,linestyle=':',color='k')
      plt.xlabel('$x$ [$\mu m$]',fontdict=font)
      plt.ylabel('$y$ [$\mu m$]',fontdict=font)
      plt.xticks([0,20,40,60],fontsize=font_size); plt.yticks([-5,0,5],fontsize=font_size);
      #plt.xscale('log')
      plt.xlim(-5, 60)
      plt.ylim(-6.5,6.5)
      plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
      plt.subplots_adjust(left=0.16, bottom=0.09, right=0.95, top=0.99, wspace=0.05, hspace=0.25)


      #plt.show()
      #lt.figure(figsize=(100,100))
      fig = plt.gcf()
      fig.set_size_inches(11, 10)
      fig.savefig(to_path+'wrap_px_py_single'+str(i).zfill(4)+'.png',format='png',dpi=160)
      plt.close("all")
      print('finished!'+str(i))
      print(np.max(gg[i,:]))
