#!/usr/bin/python

"""
calc_mem_mindcc_ncc_martini.py gro xtc
created by: Julia Rogers
created on: 10-21-20

calc min dcc, ncc for lipids (updated to include both POPC and C1P)
in membrane as function of distance from protein in xy

"""

NCCCUTOFF = 14.0 # cutoff for close contact, angstrom


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
   comtop = u.select_atoms('resid 215:346').center_of_mass(pbc=True)[2]
   combot = u.select_atoms('resid 347:478').center_of_mass(pbc=True)[2]
   comprot = prot.center_of_mass(pbc=True)[2]
   dist1 = comprot - comtop
   dist2 = combot - comprot
   box=u.dimensions
   dist1 -= round(dist1/box[2],0)*box[2]
   dist2 -= round(dist2/box[2],0)*box[2]
   if abs(dist1) < abs(dist2): topref = True
   else: topref = False
   if topref: memc = u.select_atoms('resid 215:346 and name C1A D2A C3A C4A C1B C2B C3B C4B T1A C2A C3A C1B C2B C3B C4B')
   else: memc = u.select_atoms('resid 347:478 and name C1A D2A C3A C4A C1B C2B C3B C4B T1A C2A C3A C1B C2B C3B C4B')

   print("Calculating min dcc and ncc of lipids for gro %s and xtc %s" % (sys.argv[1], sys.argv[2]))
   if topref: print("Using top leaflet as reference")
   else: print("Using bottom leaflet as reference")

   # calc min dcc, ncc as fcn of dist from prot COM
   mindccs = []
   nccs = []
   dists = []

   # number of carbons per lipid
   numc_per_lipid = np.empty(len(memc.residues), dtype=int)
   for i in range(len(memc.residues)):
      numc_per_lipid[i] = len(memc.select_atoms('resid %i' % memc.residues[i].resid))

   # loop over traj and calc quantities
   fname = 'memlip_mindcc_ncc_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for t in u.trajectory:
         if t.time % 100 != 0: continue # statistics every 100 ps
         if t.time/1000.0 % 100 == 0: print(t.time/1000.0)

         # prot center of mass
         comprot = prot.center_of_mass(pbc=True)[:2]
 
         # distances between all hydrophobic lipid carbons
         dccs = self_distance_array(memc.positions, box=t.dimensions)
         dccs2d = np.ones((len(memc),len(memc)))*1000.0
         k = 0
         for i in range(len(memc)-1):
            for j in range(i+1,len(memc)):
               dccs2d[i,j] = dccs[k]
               dccs2d[j,i] = dccs[k]
               k += 1

         # for each lipid
         for i in range(len(memc.residues)):
            comlip = u.select_atoms('resid %i' % memc.residues[i].resid).center_of_mass(pbc=True)[:2]
            # calc dist to protein
            d = calc_dist_pbc(comprot,comlip,t.dimensions[:2])
            dists.append(d)
            # get all relevant carbon-carbon dists
            sndx = np.sum(numc_per_lipid[:i])
            endx = np.sum(numc_per_lipid[:i+1])
            dcc = np.ravel(dccs2d[sndx:endx,np.r_[0:sndx,endx:len(memc)]])
            # calc min dcc and ncc
            mindcc = np.amin(dcc)
            ncc = len(np.where(dcc <= NCCCUTOFF)[0])
            mindccs.append(mindcc)
            nccs.append(ncc)
            f.write('%8.4f %5i %8.4f %8.4f %8i\n' % (t.time/1000.0,memc.residues[i].resid,d,mindcc,ncc))

   # output 2d histogram for min dcc
   hist,binx,biny = np.histogram2d(dists, mindccs, density=False, bins=[50,50])
   dxs = np.diff(binx**2)
   # normalize by dx**2 to account for area
   hist /= dxs[:,None]
   hist /= hist.sum()
   widthx = 0.5*(binx[1]-binx[0])
   widthy = 0.5*(biny[1]-biny[0])
   fname='memlip_hist2d_protdist_mindcc_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for i in range(len(hist)):
         for j in range(len(hist[0])):
            f.write('%8.4f %8.4f %16.8f\n' % (binx[i]+widthx, biny[j]+widthy, hist[i][j]))
         f.write('\n')

   # output 2d histogram for ncc
   hist,binx,biny = np.histogram2d(dists, nccs, density=False, bins=[50,50])
   dxs = np.diff(binx**2)
   # normalize by dx**2 to account for area
   hist /= dxs[:,None]
   hist /= hist.sum()
   widthx = 0.5*(binx[1]-binx[0])
   widthy = 0.5*(biny[1]-biny[0])
   fname='memlip_hist2d_protdist_ncc_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for i in range(len(hist)):
         for j in range(len(hist[0])):
            f.write('%8.4f %8.4f %16.8f\n' % (binx[i]+widthx, biny[j]+widthy, hist[i][j]))
         f.write('\n')

if __name__ == '__main__':
   main()

