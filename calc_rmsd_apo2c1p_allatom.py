#!/usr/bin/python

"""
calc_rmsd_apo2c1p_allatom.py ref_gro apo_gro apo_xtc c1p_gro c1p_xtc
created by: Julia Rogers
created on: 3-20-22

calculates av CA rmsd between apo and c1p structures

"""


import sys
import numpy as np
from scipy.spatial.distance import cdist
from MDAnalysis import *
from MDAnalysis.analysis.rms import RMSD, rmsd
from MDAnalysis.analysis import align
from MDAnalysis.analysis.distances import dist
from helper.general import get_basename

nres = 214

def main():
   uapo = Universe(sys.argv[2], sys.argv[3])
   uc1p = Universe(sys.argv[4], sys.argv[5])

   prot = uapo.select_atoms('protein')
   if prot.resids[0] > 1: md_resid_offset = 264
   else: md_resid_offset = 0

   # align based on helix 6
   selection = 'resid %i-%i and name CA' % (152+md_resid_offset, 163+md_resid_offset)
   # use frame 1 from c1p as reference
   ref = Universe(sys.argv[1])
   # align c1p traj
   aligner = align.AlignTraj(uc1p, ref, select=selection, in_memory=True).run()
   # align apo traj
   aligner = align.AlignTraj(uapo, ref, select=selection, in_memory=True).run() 

   carmsds = np.zeros(nres)
   with open('apo_c1p_carmds.txt', 'w') as f:
      for i in range(nres):
         selection = 'resid %i and name CA' % (i+1+md_resid_offset)
         aposel = uapo.select_atoms(selection)
         c1psel = uc1p.select_atoms(selection)
         apocoords = np.array([aposel.positions[0] for t in uapo.trajectory])
         c1pcoords = np.array([c1psel.positions[0] for t in uc1p.trajectory])
         dists = cdist(apocoords, c1pcoords, metric='euclidean')
         dists = np.power(dists, 2)
         carmsds[i] = np.sum(dists)/len(uapo.trajectory)/len(uc1p.trajectory)
         carmsds[i] = np.sqrt(carmsds[i])
         f.write('%5i %8.4f\n' % ((i+1), carmsds[i]))


if __name__ == '__main__':
   main()

