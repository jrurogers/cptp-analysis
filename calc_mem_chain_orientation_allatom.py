#!/usr/bin/python

"""
calc_mem_chain_orientation_allatom.py gro xtc
created by: Julia Rogers
created on: 12-1-21

calc cos theta_z of vector to each lipid chain
oriented towards the center of the bilayer: cos theta_z = 1

"""


import sys, os, re
import numpy as np
from MDAnalysis.analysis.distances import distance_array, self_distance_array
from MDAnalysis.analysis.contacts import soft_cut_q
from MDAnalysis.lib.util import convert_aa_code
from MDAnalysis import *

def calc_dist_pbc(a,b,box):
   d = a - b
   d -= np.around(d/box)*box
   return np.sqrt(np.sum(d*d))

def main():

   # select atoms
   u=Universe(sys.argv[1], sys.argv[2])
   prot = u.select_atoms('protein')

   # determine ref leaflet that protein is bound to
   comtop = u.select_atoms('resid 1:132').center_of_mass(pbc=True)[2]
   combot = u.select_atoms('resid 133:264').center_of_mass(pbc=True)[2]
   comprot = prot.center_of_mass(pbc=True)[2]
   dist1 = comprot - comtop
   dist2 = combot - comprot
   box=u.dimensions
   dist1 -= round(dist1/box[2],0)*box[2]
   dist2 -= round(dist2/box[2],0)*box[2]
   if abs(dist1) < abs(dist2): topref = True
   else: topref = False

   if topref: resid_sel='1:132'
   else: resid_sel='133:264'
   memc2 = u.select_atoms('resid %s and name C2 C2S' % resid_sel)
   memt = u.select_atoms('resid %s and name C218 C316 C18S C16F' % resid_sel)

   print "Calculating orientation of lipid tails for gro %s and xtc %s" % (sys.argv[1], sys.argv[2])
   if topref: print "Using top leaflet as reference"
   else: print "Using bottom leaflet as reference"

   # calc cos theta of tail vector with z axis
   costhetas = []
   dists = []

   # loop over traj and calc quantities
   fname = 'memlip_costhetaz_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for t in u.trajectory:
         if t.time/1000.0 < 50: continue # statistics from 50-100 ns
         if t.time/1000.0 % 10 == 0: print t.time/1000.0

         # prot center of mass
         comprot = prot.center_of_mass(pbc=True)[:2]

         # tail vectors
         if topref:
            vec1 = memc2.positions - memt.positions[::2]
            vec2 = memc2.positions - memt.positions[1::2]
         else:
            vec1 = memt.positions[::2] - memc2.positions
            vec2 = memt.positions[1::2] - memc2.positions

         # normalized tail vectors
         vec1 /= np.sqrt(np.einsum('...i,...i', vec1, vec1))[...,np.newaxis]
         vec2 /= np.sqrt(np.einsum('...i,...i', vec2, vec2))[...,np.newaxis]
         
         # for each lipid
         for i in range(len(memc2.residues)):
            comlip = u.select_atoms('resid %i' % memc2.residues[i].resid).center_of_mass(pbc=True)[:2]
            # calc dist to protein
            d = calc_dist_pbc(comprot,comlip,t.dimensions[:2])
            dists.append(d)
            dists.append(d)
            costhetas.append(vec1[i,2])
            costhetas.append(vec2[i,2])
            
            f.write('%8.4f %5i %8.4f %8.4f %8.4f\n' % (t.time/1000.0,memc2.residues[i].resid,d,vec1[i,2],vec2[i,2]))

   # output 2d histogram for cos theta
   hist,binx,biny = np.histogram2d(dists, costhetas, density=False, bins=[50,50])
   dxs = np.diff(binx**2)
   # normalize by dx**2 to account for area
   hist /= dxs[:,None]
   hist /= hist.sum()
   widthx = 0.5*(binx[1]-binx[0])
   widthy = 0.5*(biny[1]-biny[0])
   fname='memlip_hist2d_protdist_costhetaz_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for i in range(len(hist)):
         for j in range(len(hist[0])):
            f.write('%8.4f %8.4f %16.8f\n' % (binx[i]+widthx, biny[j]+widthy, hist[i][j]))
         f.write('\n')


if __name__ == '__main__':
   main()

