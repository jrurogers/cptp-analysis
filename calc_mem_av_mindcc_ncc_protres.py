#!/usr/bin/python

"""
calc_mem_av_mindcc_ncc_protres.py gro xtc memlip_mindcc_ncc.txt
created by: Julia Rogers
created on: 7-16-20

calc av min dcc and ncc of lipids w/heavy atom w/in 5 ang
of protein residue heavy atom


* note: must have already calculated min dcc and ncc for each lipid *

"""


import sys, os
import numpy as np
from MDAnalysis import *
from MDAnalysis import transformations
from MDAnalysis.analysis.distances import distance_array, self_distance_array
from MDAnalysis.analysis.density import density_from_Universe

def calc_dist_pbc(a,b,box):
   d = a - b
   d -= np.around(d/box)*box
   return np.sqrt(np.sum(d*d))

def wrap(pos,box):
   pos -= np.floor(pos/box)*box
   return pos

def main():
   u = Universe(sys.argv[1], sys.argv[2])
   prot = u.select_atoms('protein')
   # determine ref leaflet that protein is bound to
   comtop = u.select_atoms('resid 1:132').center_of_mass(pbc=True)[2]
   combot = u.select_atoms('resid 133:264').center_of_mass(pbc=True)[2]
   comprot = prot.center_of_mass(pbc=True)
   box=u.dimensions
   dist1 = comprot[2] - comtop
   dist2 = combot - comprot[2]
   dist1 -= round(dist1/box[2],0)*box[2]
   dist2 -= round(dist2/box[2],0)*box[2]
   if abs(dist1) < abs(dist2): topref = True
   else: topref = False
   if topref: memsel = 'resid 1:132 and not name H*'
   else: memsel = 'resid 133:264 and not name H*'

   print "Calculating av min dcc and ncc of lipids for each prot residue for gro %s and xtc %s" % (sys.argv[1], sys.argv[2])
   if topref: print "Using top leaflet as reference"
   else: print "Using bottom leaflet as reference"

   # protein
   prot = u.select_atoms('protein and not name H*')
   nres = len(prot.residues)
   print(nres)
   # number of atoms per residue
   numa_per_residue = np.empty(nres, dtype=int)
   for i in range(nres):
      numa_per_residue[i] = len(prot.select_atoms('resid %i' % prot.residues[i].resid))

   # lipids
   mem = u.select_atoms(memsel)
   nlip = len(mem.residues)
   # number of atoms per lipid
   numa_per_lipid = np.empty(nlip, dtype=int)
   for i in range(nlip):
      numa_per_lipid[i] = len(mem.select_atoms('resid %i' % mem.residues[i].resid))

   # get min dcc and ncc data
   mindccs = np.reshape(np.loadtxt(sys.argv[3], usecols=(3,)), (len(u.trajectory), nlip))
   nccs = np.reshape(np.loadtxt(sys.argv[3], usecols=(4,)), (len(u.trajectory), nlip))
   print(mindccs.shape, nccs.shape)

   res_mindccs = [[] for i in range(nres)]
   res_nccs = [[] for i in range(nres)]
   fname = 'protres_nearby_lipid_mindcc_ncc_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for k, t in enumerate(u.trajectory):
         if t.time/1000.0 < 50: continue # statistics from 50-100 ns
         if t.time/1000.0 % 10 == 0: print t.time/1000.0

         # distances btwn all heavy atoms of lipids and prot residues
         ds = distance_array(mem.positions, prot.positions, box=t.dimensions)

         for i in range(nres):
            curr_mindcc = []
            curr_ncc = []
            # get all relevant prot res dists
            sndx = np.sum(numa_per_residue[:i])
            endx = np.sum(numa_per_residue[:i+1])
            # see if each lipid has min dist <= 5 ang
            for j in range(nlip):
               lsndx = np.sum(numa_per_lipid[:j])
               lendx = np.sum(numa_per_lipid[:j+1])
               mind = np.amin(ds[lsndx:lendx,sndx:endx])
               if mind <= 5:
                  curr_mindcc.append(mindccs[k,j])
                  curr_ncc.append(nccs[k,j])
            res_mindccs[i].extend(curr_mindcc)
            res_nccs[i].extend(curr_ncc)
            f.write('%8.4f %5i %8.4f %8.4f\n' % (t.time/1000.0, i+1, np.mean(curr_mindcc), np.mean(curr_ncc)))


   fname = 'memlip_av_mindcc_ncc_protres_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for i in range(nres):
         if len(res_nccs[i]) == 0: f.write('%5i -1 -1 -1 -1\n' % (i+1))
         else: f.write('%5i %8.4f %8.4f %8.4f %8.4f\n' % (i+1, np.mean(res_mindccs[i]), np.mean(res_nccs[i]), np.std(res_mindccs[i]), np.std(res_nccs[i])))


if __name__ == '__main__':
   main()

