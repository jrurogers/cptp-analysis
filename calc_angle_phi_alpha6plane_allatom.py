#!/usr/bin/python

"""
calc_angle_phi_alpha6plane_allatom.py gro xtc
created by: Julia Rogers
created on: 11-13-20

calculates angle between helix 2 and normal of plane of helix 6

"""

NBINS=90

import sys
import numpy as np
from helper.general import get_basename
from math import degrees,atan2
from MDAnalysis import *
from MDAnalysis.lib.mdamath import angle, dihedral

# returns unit vector normal to plane defining surface of protein
def define_plane(c,n,ca):
   z = np.cross((ca-c), (n-c))
   z = z / np.linalg.norm(z)
   return z

def main():
   u=Universe(sys.argv[1], sys.argv[2])
   prot = u.select_atoms('protein')
   if prot.resids[0] > 1: memoffset = 264
   else: memoffset = 0

   outname = get_basename(sys.argv[2])[:-4]

   # atoms to define plane of helix 6
   ca155 = u.select_atoms('resid %i and name CA' % (155+memoffset))
   ca158 = u.select_atoms('resid %i and name CA' % (158+memoffset))
   ca162 = u.select_atoms('resid %i and name CA' % (162+memoffset))

   # atoms to define vector of helix 2
   ca54 = u.select_atoms('resid %i and name CA' % (54+memoffset))
   ca65 = u.select_atoms('resid %i and name CA' % (65+memoffset))

   avtheta=0.0
   avphi=0.0

   thetas=[]
   phis=[]

   with open("angle_phi_alpha6plane_"+outname+".txt", 'w') as f:
      for t in u.trajectory:
         if t.time/1000.0 % 10 == 0: print t.time/1000.0
         # coordinates of points to define plane
         a = np.ravel(ca155.positions)
         b = np.ravel(ca158.positions)
         c = np.ravel(ca162.positions)
         # determine vector normal to plane
         z = define_plane(a,b,c)

         # helix 2 vector
         vec = np.ravel(ca65.positions) - np.ravel(ca54.positions)
         # angle between vector normal to plane and normal vector
         theta = degrees(angle(vec,z))
         # angle between helix 6 axis and projection of helix 2 onto plane
         helix = a - c
         proj = vec - np.dot(vec,z) * z
         dot = np.dot(proj, helix)
         det = np.dot(z, np.cross(proj, helix))
         phi = degrees(atan2(det, dot))
     
         # calc average angles
         avtheta += theta
         avphi += phi

         # output angles at each time
         thetas.append(theta)
         phis.append(phi)
         f.write('%8.4f %8.4f %8.4f\n' % (t.time/1000.0, theta, phi))

   avtheta /= float(len(u.trajectory))
   avphi /= float(len(u.trajectory))
      
   with open("av_angle_phi_alpha6plane_"+outname+".txt", 'w') as f:
      f.write("average angle: %8.4f\n" % avtheta)
      f.write("average phi: %8.4f\n" % avphi)

   hist, bins = np.histogram(thetas, range=[0,180], bins=NBINS, density=True)
   width=(bins[1]-bins[0])/2.0
   with open("distrib_angle_alpha6plane_"+outname+".txt", 'w') as f:
      for i in range(len(hist)):
         f.write('%8.4f %10.6f\n' % ((bins[i]+width), hist[i]))

   hist, bins = np.histogram(phis, range=[-180,180], bins=NBINS, density=True)
   width=(bins[1]-bins[0])/2.0
   with open("distrib_phi_alpha6plane_"+outname+".txt", 'w') as f:
      for i in range(len(hist)):
         f.write('%8.4f %10.6f\n' % ((bins[i]+width), hist[i]))

   hist, xbins, ybins = np.histogram2d(thetas, phis, range=[[0,180], [-180,180]], bins=[90,90], density=True)
   xwidth = (xbins[1]-xbins[0])/2.0
   ywidth = (ybins[1]-ybins[0])/2.0
   with open("distrib2d_theta_phi_refcoord_"+outname+".txt", 'w') as f:
      for i in range(len(hist)):
         for j in range(len(hist[i])):
            f.write('%8.4f %8.4f %12.8f\n' % (xbins[i]+xwidth, ybins[j]+ywidth, hist[i,j]))
         f.write('\n')




if __name__ == '__main__':
   main()

