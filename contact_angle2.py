
#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
http://www.scipy.org/Cookbook/Least_Squares_Circle
"""

from numpy import *


import numpy as np
import matplotlib.pyplot as plt
import sys
from dump import dump
from scipy.optimize import curve_fit
import math
from math import sqrt, atan, pi


d=dump("dump3.lammpstrj")

d.map(1,"id",2,"type",3,"x",4,"y",5,"z")

d.tselect.test("$t >=16920000 and $t <=20000000")  # the frames of the lammps Trajectory file
d.delete()
d.sort()
d.aselect.test("$type!=14") # The atom type of the substrate
d.write("temp/tmp.lammpstrj")
d = dump("temp/tmp.lammpstrj",0)
d.map(1,"id",2,"type",3,"x",4,"y",5,"z")


#time = d.next()
#print 'time', time
rbin = 100   #number of radial bins
zbin = 100  #number of z bins
atom = np.zeros((zbin,rbin)) 
xdata = np.zeros((zbin-1))
z4 = []
r = [] #np.zeros((zbin-1))
ydata = np.zeros((rbin-1))
zdata = np.zeros((rbin-1))

da = 95
nsnaps = 0
success = True
outfile = "temp/out.dat"
outfile2 = "temp/out2.dat"
while 1:
  time = d.next()
  #d.scale(time)
  if time == -1: break
  #print time

  x1 = d.vecs(time,"x")
  y1 = d.vecs(time,"y")
  z1 = d.vecs(time,"z")
  natom = len(x1)
  minz,maxz = d.minmax("z")
  zl = minz + 3.35*2   #5  #min(z1) + 20 #+ 3.4    # +3.4*n (added for pillar height, n = height of pillar) # The reference height 
  zu = minz + zbin
  xcom = sum(x1)/natom
  ycom = sum(y1)/natom
  shiftz = 5
  #dist = len(x)
  #dist = ((x-xcom)**2 + (y-ycom)**2)**(1/2)
  for i in xrange(natom):
      #dist = ((x1[i]-xcom)**2 #+ (y1[i]-ycom)**2)
      #dist = dist**0.5
      #print 'x and xcom and y and ycom', i, x[i], xcom, y[i], ycom
      #dist1 = (rbin*da/3.14)
      #dist1 = dist1**0.5
      #print 'dist and dist1', dist, dist1
      #if(dist <= dist1):
      #  j = 0
      #  while success:
      #    j = j + 1
      #    dist2= ((j*da/3.14))**0.5
      #    if(dist <= dist2): success = False     
          #print 'dist and dist2', dist, dist2           
      #success = True
      x2 = abs(y1[i] - ycom)
      x2 = int(x2) + 1
      z2 = z1[i] - zl
      z2 = int(z2) + 1
      z2 = z2 -shiftz
      #if(z1 > 200): print ' z and z1', i, z[i], z1
      #print 'z1 j', z1, j
      if(z2 > 0 and z2 < zbin-1 and x2 < rbin-1): atom[z2,x2] = atom[z2,x2] + 1



  nsnaps += 1
#  print time
atom = atom/nsnaps
#print 'atom', atom
#*****************************************printing radial and axial density*****************

fp = open(outfile,"w")
for i in xrange(zbin-1):
  for j in xrange(rbin-1):
    print >>fp, zl+i+0.5, xcom+j+0.5, atom[i+1,j+1],atom[i+1,j+1]/3.1
  print >>fp
fp.close()

#***************************************finding the interface points******************************

def func(x, a, b, d, r):
    return 0.5*(a+b) - 0.5*(a-b)*np.tanh((2*(x-r)/d))

# Summary
fmt = '%-22s %20.5f %20.5f %20.6f %20.6f'
print ('\n%-22s' +' %20s'*4) % tuple('METHOD liq-density vapor-denisity interfacial-width interfacial-points'.split())
print '-'*(22 +7*(10+1))
method  = "tan hyperbolic"

fp = open(outfile2,"w")
from matplotlib                 import pyplot as p, cm, colors
f = p.figure( facecolor='white')  #figsize=(7, 5.4), dpi=72,
p.close('all')
for i in xrange(zbin-1):
  xdata[i] = i+0.5+shiftz
  for j in xrange(rbin-1):
    ydata[j] = j + 0.5
    zdata[j] = atom[i+1,j+1]
  #print 'ydata and zdata', ydata,zdata 
  popt, pcov = curve_fit(func, ydata, zdata)
  print >> fp, popt
  r.extend([popt[3]])
  z4.extend([xdata[i]])
  if popt[3]<= 12: break
  print  fmt % (method, popt[0], popt[1], popt[2], popt[3] )
  if i == 5: 
    zfit = func(ydata,*popt)
    p.plot(ydata,zdata, '-', label='method_tanh', lw=2)
    p.plot(ydata,zfit, '*', label='method_tanh', lw=2)
    p.savefig('fig1.png')
    p.xlabel('Radial coordinate')
    p.ylabel('Density')
    p.close(f)
  
fp.close()



#************************************fitting interface****************************************************
#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
http://www.scipy.org/Cookbook/Least_Squares_Circle
"""

from numpy import *


y = z4 
x = r
x_m = mean(x)
y_m = mean(y)



# Decorator to count functions calls
import functools
def countcalls(fn):
    "decorator function count function calls "

    @functools.wraps(fn)
    def wrapped(*args):
        wrapped.ncalls +=1
        return fn(*args)

    wrapped.ncalls = 0
    return wrapped



from scipy      import optimize

method_2  = "leastsq"

def calc_R(xc, yc):
    """
        # calculate the distance of each 2D points from the center (xc, yc) 
    """
    return sqrt((x-xc)**2 + (y-yc)**2)

@countcalls
def f_2(c):
    """ 
    # calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) 
    """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

center_estimate = x_m, y_m
center_2, ier = optimize.leastsq(f_2, center_estimate)

xc_2, yc_2 = center_2
Ri_2       = calc_R(xc_2, yc_2)
R_2        = Ri_2.mean()
residu_2   = sum((Ri_2 - R_2)**2)
residu2_2  = sum((Ri_2**2-R_2**2)**2)
ncalls_2   = f_2.ncalls
basename = 'circle'
print 'center', center_2, R_2
x0= sqrt(R_2**2 - yc_2**2)
slope = -x0/yc_2
angle = (180*math.atan(slope))/pi
print 'angle', angle



def plot_all(residu2=False):
    """
    #Draw data points, best fit circles and center for the three methods,
    #and adds the iso contours corresponding to the fiel residu or residu2
    """

    f = p.figure( facecolor='white')  #figsize=(7, 5.4), dpi=72,
    p.axis('equal')

    theta_fit = linspace(-pi, pi, 180)


    x_fit2 = xc_2 + R_2*cos(theta_fit)
    y_fit2 = yc_2 + R_2*sin(theta_fit)
    p.plot(x_fit2, y_fit2, 'k--', label=method_2, lw=2)


    p.plot([xc_2], [yc_2], 'gD', mec='r', mew=1)

    # draw
    p.xlabel('x')
    p.ylabel('y')

    # plot the residu fields
    nb_pts = 100

    p.draw()
    xmin, xmax = p.xlim()
    ymin, ymax = p.ylim()

    vmin = min(xmin, ymin)
    vmax = max(xmax, ymax)

    xg, yg = ogrid[vmin:vmax:nb_pts*1j, vmin:vmax:nb_pts*1j]
    xg = xg[..., newaxis]
    yg = yg[..., newaxis]

    Rig    = sqrt( (xg - x)**2 + (yg - y)**2 )
    Rig_m  = Rig.mean(axis=2)[..., newaxis]

    residu = sum( (Rig**2 - Rig_m**2)**2 ,axis=2)

    lvl = exp(linspace(log(residu.min()), log(residu.max()), 15))

    p.contourf(xg.flat, yg.flat, residu.T, lvl, alpha=0.4, cmap=cm.Purples_r) # , norm=colors.LogNorm())
    cbar = p.colorbar(fraction=0.175, format='%.f')
    p.contour (xg.flat, yg.flat, residu.T, lvl, alpha=0.8, colors="lightblue")

    cbar.set_label('Residu_2 - algebraic approximation')

    # plot data
    p.plot(x, y, 'ro', label='data', ms=8, mec='b', mew=1)
    p.legend(loc='best',labelspacing=0.1 )

    p.xlim(xmin=vmin, xmax=vmax)
    p.ylim(ymin=vmin, ymax=vmax)

    p.grid()
    p.title('Least Squares Circle')
    p.savefig('%s_residu%d.png' % (basename, 2 if residu2 else 1))

#plot_all(residu2=False)
plot_all(residu2=True )

p.show()

