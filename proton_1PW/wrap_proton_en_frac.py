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
        'size'   : 20,  
        }  
 
font_size=25
font_size2=20

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
  part_name='wrap_bar_proton' 

  plt.subplot(1,3,1)
  from_path = './PW_w010/'
  to_path   = './PW_w010_fig/'
  px_d = np.loadtxt(from_path+part_name+'_px.txt')
  py_d = np.loadtxt(from_path+part_name+'_py.txt')
  xx_d = np.loadtxt(from_path+part_name+'_xx.txt')
  yy_d = np.loadtxt(from_path+part_name+'_yy.txt')
  ek_d = ((px_d**2+py_d**2+1)**0.5-1.0)*0.51*mass  
  value_axisx = np.linspace(0,340,34)
  value_grid = np.linspace(0,340,35)
  value_left = np.zeros_like(value_axisx)
  value_mid  = np.zeros_like(value_axisx)
  value_right= np.zeros_like(value_axisx)
  for i in range(np.size(value_axisx)):
        final_ek       = ek_d[:,1]
        initial_x      = xx_d[:,0]
        where_in_bin   = np.where((value_grid[i]<=final_ek) & (value_grid[i+1]>final_ek))
        print(where_in_bin)
        in_x           = initial_x[where_in_bin]
        value_left[i]  = np.size(in_x[(in_x<=0.1) & (in_x>0)])
        value_mid[i]   = np.size(in_x[(in_x<=0.2) & (in_x>0.1)])
        value_right[i] = np.size(in_x[(in_x<=0.3) & (in_x>0.2)])
  y_l = value_left/(value_left+value_mid+value_right)
  y_m = value_mid/(value_left+value_mid+value_right)
  y_r = value_right/(value_left+value_mid+value_right)
  print(y_l,y_m,y_r)
  width=10
  pl=plt.bar(value_axisx, y_l*100, width, color='deepskyblue',edgecolor='black',linewidth=1)
  pt=plt.bar(value_axisx, y_m*100, width, bottom=y_l*100, color='lightgreen',edgecolor='black',linewidth=1)
  pt=plt.bar(value_axisx, y_r*100, width, bottom=(y_l+y_m)*100, color='salmon',edgecolor='black',linewidth=1)
  plt.xlim(-5,335)
  plt.ylim(0,100)
  plt.xlabel('$\epsilon_p$ [MeV]',fontdict=font)
  plt.ylabel('N$_{l,m,r}$/N$_{tot}$ [%]',fontdict=font)
  plt.xticks(fontsize=font_size);
  plt.yticks(fontsize=font_size);
#  plt.title('proton initial fraction',fontdict=font)


  plt.subplot(1,3,2)
  from_path = './PW_w020/'
  to_path   = './PW_w020_fig/'
  px_d = np.loadtxt(from_path+part_name+'_px.txt')
  py_d = np.loadtxt(from_path+part_name+'_py.txt')
  xx_d = np.loadtxt(from_path+part_name+'_xx.txt')
  yy_d = np.loadtxt(from_path+part_name+'_yy.txt')
  ek_d = ((px_d**2+py_d**2+1)**0.5-1.0)*0.51*mass  
  value_axisx = np.linspace(0,230,34)
  value_grid = np.linspace(0,230,35)
  value_left = np.zeros_like(value_axisx)
  value_mid  = np.zeros_like(value_axisx)
  value_right= np.zeros_like(value_axisx)
  for i in range(np.size(value_axisx)):
        final_ek       = ek_d[:,1]
        initial_x      = xx_d[:,0]
        where_in_bin   = np.where((value_grid[i]<=final_ek) & (value_grid[i+1]>final_ek))
        print(where_in_bin)
        in_x           = initial_x[where_in_bin]
        value_left[i]  = np.size(in_x[(in_x<=0.1) & (in_x>0)])
        value_mid[i]   = np.size(in_x[(in_x<=0.2) & (in_x>0.1)])
        value_right[i] = np.size(in_x[(in_x<=0.3) & (in_x>0.2)])
  y_l = value_left/(value_left+value_mid+value_right)
  y_m = value_mid/(value_left+value_mid+value_right)
  y_r = value_right/(value_left+value_mid+value_right)
  print(y_l,y_m,y_r)
  width=230./34
  pl=plt.bar(value_axisx, y_l*100, width, color='deepskyblue',edgecolor='black',linewidth=1)
  pt=plt.bar(value_axisx, y_m*100, width, bottom=y_l*100, color='lightgreen',edgecolor='black',linewidth=1)
  pt=plt.bar(value_axisx, y_r*100, width, bottom=(y_l+y_m)*100, color='salmon',edgecolor='black',linewidth=1)
  plt.xlim(-width/2,225)
  plt.ylim(0,100)
  plt.xlabel('$\epsilon_p$ [MeV]',fontdict=font)
#  plt.ylabel('N$_{l,m,r}$/N$_{tot}$ [%]',fontdict=font)
  plt.xticks(fontsize=font_size);
  plt.yticks(fontsize=0.001);
#  plt.title('proton initial fraction',fontdict=font)



  plt.subplot(1,3,3)
  from_path = './PW_w080/'
  to_path   = './PW_w080_fig/'
  px_d = np.loadtxt(from_path+part_name+'_px.txt')
  py_d = np.loadtxt(from_path+part_name+'_py.txt')
  xx_d = np.loadtxt(from_path+part_name+'_xx.txt')
  yy_d = np.loadtxt(from_path+part_name+'_yy.txt')
  ek_d = ((px_d**2+py_d**2+1)**0.5-1.0)*0.51*mass  
  value_axisx = np.linspace(0,50,34)
  value_grid = np.linspace(0,50,35)
  value_left = np.zeros_like(value_axisx)
  value_mid  = np.zeros_like(value_axisx)
  value_right= np.zeros_like(value_axisx)
  for i in range(np.size(value_axisx)):
        final_ek       = ek_d[:,1]
        initial_x      = xx_d[:,0]
        where_in_bin   = np.where((value_grid[i]<=final_ek) & (value_grid[i+1]>final_ek))
        print(where_in_bin)
        in_x           = initial_x[where_in_bin]
        value_left[i]  = np.size(in_x[(in_x<=0.1) & (in_x>0)])
        value_mid[i]   = np.size(in_x[(in_x<=0.2) & (in_x>0.1)])
        value_right[i] = np.size(in_x[(in_x<=0.3) & (in_x>0.2)])
  y_l = value_left/(value_left+value_mid+value_right)
  y_m = value_mid/(value_left+value_mid+value_right)
  y_r = value_right/(value_left+value_mid+value_right)
  print(y_l,y_m,y_r)
  width=50/34
  pl=plt.bar(value_axisx, y_l*100, width, color='deepskyblue',edgecolor='black',linewidth=1)
  pt=plt.bar(value_axisx, y_m*100, width, bottom=y_l*100, color='lightgreen',edgecolor='black',linewidth=1)
  pt=plt.bar(value_axisx, y_r*100, width, bottom=(y_l+y_m)*100, color='salmon',edgecolor='black',linewidth=1)
  plt.xlim(-width/2,48)
  plt.ylim(0,100)
  plt.xlabel('$\epsilon_p$ [MeV]',fontdict=font)
#  plt.ylabel('N$_{l,m,r}$/N$_{tot}$ [%]',fontdict=font)
  plt.xticks(fontsize=font_size);
  plt.yticks(fontsize=0.001);
#  plt.legend(['left','mid','right'],loc='best',fontsize=font_size,framealpha=0.5)
#  plt.title('proton initial fraction',fontdict=font)



#  plt.grid(which='major',linestyle=':', linewidth=1, color='k')

  plt.subplots_adjust(left=0.14, bottom=0.22, right=0.98, top=0.95, wspace=0.08, hspace=0.01)
  fig = plt.gcf()
  fig.set_size_inches(10, 4.)
  fig.savefig('wrap_proton_en_frac.png',format='png',dpi=160)
  plt.close("all")
  print('wrap_proton_en_frac.png')

