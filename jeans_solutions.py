from scipy import interpolate 
import scipy
#from pydl import smooth
import numpy as np
import os.path
import sys

#it works only with scipy version >= 0.17.0

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def jeans_equation_solution():
  print("\nVelocites obtained solving Jean's Equations for multiple stellar population\n")
  g = open("inputJE.txt","r")
  if os.path.exists("inputJE.txt") == False:
    sys.exit("File inputJE.txt not found")
  else:
    print("inputJE.txt opened successfully")


  ntot = int(g.readline().rstrip('\n'))
  seed = int(g.readline().rstrip('\n'))

  lines=g.readlines()
  rho = []
  m = []
  r = []
  for x in lines:
    r.append(float(x.split(' ')[0]))
    rho.append(float(x.split(' ')[1]))
    m.append(float(x.split(' ')[2]))
  g.close()
  rho = np.array(rho); r = np.array(r); m = np.array(m)
  ngrd=30000
  sigma=np.zeros(ngrd)
  aa=np.arange(float(ngrd+1))/ngrd
  rmin=min(r)
  rmax=max(r)
  rgrd=np.exp((np.log(rmax)-np.log(rmin))*aa+np.log(rmin))
  rho=smooth(rho,10) #Rho is the density profiles
  interpfuc = interpolate.interp1d(r,rho,kind='linear',fill_value='extrapolate')
  rhogrd = interpfuc(rgrd)
#  for i in range(0,ngrd):
#    if rhogrd[i] == 0.0:
#      rhogrd[i] = 1E-8 # taking care of the no tidal field case: rho could have been extrapolated with zero values

  interpfuc=interpolate.interp1d(r,m,kind='linear',fill_value='extrapolate') #m is the cumulative mass profiles
  mgrd = interpfuc(rgrd)

  dwdr=mgrd/rgrd**2
  wgrd=np.zeros(ngrd)
  for i in range(0,ngrd):
    wgrd[i]=np.trapz(dwdr[i:],rgrd[i:])
    
  #we can obtain rho and m within mcluster with the function radial_profile
  for i in range(0,ngrd):
    sigma[i]=1.0/rhogrd[i]*np.trapz(rhogrd[i:]*dwdr[i:],rgrd[i:])

  interpfuc=interpolate.interp1d(rgrd[0:ngrd],np.sqrt(sigma),kind='linear',fill_value='extrapolate')
  veld = interpfuc(r)

  interpfuc=interpolate.interp1d(rgrd[0:ngrd],np.sqrt(2*wgrd),kind='linear',fill_value='extrapolate')
  vesc = interpfuc(r)

  vx=np.zeros(ntot)
  vy=np.zeros(ntot)
  vz=np.zeros(ntot)
  np.random.seed(seed)
  for i in range(0,ntot):
    V = vesc[i]
    while V > 0.99*vesc[i]:
      vx[i]=veld[i]*np.random.standard_normal()
      vy[i]=veld[i]*np.random.standard_normal()
      vz[i]=veld[i]*np.random.standard_normal()
      V = np.sqrt(vx[i]**2+vy[i]**2+vz[i]**2)
    #print(vx[i],vy[i],vz[i])
#  print(vx)
  f = open('outputJE.txt',"w")
  for i in range(0,ntot):
    f.write("{0:.6e} {1:.6e} {2:.6e}\n".format(vx[i],vy[i],vz[i]))
  f.close()

if __name__ == "__main__":
 jeans_equation_solution ()
