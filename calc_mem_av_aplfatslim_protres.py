#!/usr/bin/python

"""
calc_mem_av_aplfatslim_protres.py gro xtc apl_basename
created by: Julia Rogers
created on: 11-9-22

calc av APL of lipids w/heavy atom w/in 5 ang
of protein residue heavy atom

* note: must have already calculated APL for each lipid *

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

def load_fatslim_apl(fbase, resids, nframes):
   nres = len(resids)
   apls = np.zeros((nframes, nres))

   for i in range(nframes):
      d = np.loadtxt('%s_frame_%05d.csv' % (fbase, i), delimiter=',', usecols=(0,5), skiprows=1)
      mask = np.isin(d[:,0], resids)
      apls[i,:] = d[mask,1]*100.0
   return apls

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

   print "Calculating av APL of lipids for each prot residue for gro %s and xtc %s" % (sys.argv[1], sys.argv[2])
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

   # get apl data
   if topref: resids = [i for i in range(1,133)]
   else: resids = [i for i in range(133,265)]
   apls = load_fatslim_apl(sys.argv[3], resids, len(u.trajectory))
   print(apls.shape)
   print(apls[0,0])

   res_vals = [[] for i in range(nres)]
   fname = 'protres_nearby_lipid_aplfatslim_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for k, t in enumerate(u.trajectory):
         if t.time/1000.0 < 50: continue # statistics from 50-100 ns
         if t.time/1000.0 % 10 == 0: print t.time/1000.0

         # distances btwn all heavy atoms of lipids and prot residues
         ds = distance_array(mem.positions, prot.positions, box=t.dimensions)

         for i in range(nres):
            curr = []
            # get all relevant prot res dists
            sndx = np.sum(numa_per_residue[:i])
            endx = np.sum(numa_per_residue[:i+1])
            # see if each lipid has min dist <= 5 ang
            for j in range(nlip):
               lsndx = np.sum(numa_per_lipid[:j])
               lendx = np.sum(numa_per_lipid[:j+1])
               mind = np.amin(ds[lsndx:lendx,sndx:endx])
               if mind <= 5:
                  curr.append(apls[k,j])
            res_vals[i].extend(curr)
            f.write('%8.4f %5i %8.4f\n' % (t.time/1000.0, i+1, np.mean(curr)))


   fname = 'memlip_av_aplfatslim_protres_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for i in range(nres):
         if len(res_vals[i]) == 0: f.write('%5i -1 -1\n' % (i+1))
         else: f.write('%5i %8.4f %10.6f\n' % (i+1, np.mean(res_vals[i]), np.std(res_vals[i])))


if __name__ == '__main__':
   main()

