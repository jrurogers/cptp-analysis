#!/usr/bin/python

"""
calc_rmsd2pdb_martini.py gro xtc pdb
created by: Julia Rogers
created on: 11-4-19

calculates rmsd between martini protein backbone and pdb CA reference
also calculates for each helix individually

"""


import sys
import numpy as np
from MDAnalysis import *
from MDAnalysis.analysis.rms import RMSD, rmsd
from MDAnalysis.analysis import align
from helper.general import get_basename


def main():
   u=Universe(sys.argv[1], sys.argv[2])
   sel=u.select_atoms('(resid 8-214) and name BB')
   ref=Universe(sys.argv[3])
   refsel=ref.select_atoms('(resid 8-214) and name BB')

   # helix definitions
   helices = {'alphaN': u.select_atoms('resid 10-20 and name BB'),\
              'alpha1': u.select_atoms('resid 29-46 and name BB'),\
              'alpha2': u.select_atoms('resid 51-69 and name BB'),\
              'alpha3': u.select_atoms('resid 79-89 and name BB'),\
              'alpha4': u.select_atoms('resid 104-127 and name BB'),\
              'alpha5': u.select_atoms('resid 134-144 and name BB'),\
              'alpha6': u.select_atoms('resid 152-163 and name BB'),\
              'alpha7': u.select_atoms('resid 168-174 and name BB'),\
              'alpha8': u.select_atoms('resid 180-208 and name BB')}
   refhelices = {'alphaN': ref.select_atoms('resid 10-20 and name BB'),\
                 'alpha1': ref.select_atoms('resid 29-46 and name BB'),\
                 'alpha2': ref.select_atoms('resid 51-69 and name BB'),\
                 'alpha3': ref.select_atoms('resid 79-89 and name BB'),\
                 'alpha4': ref.select_atoms('resid 104-127 and name BB'),\
                 'alpha5': ref.select_atoms('resid 134-144 and name BB'),\
                 'alpha6': ref.select_atoms('resid 152-163 and name BB'),\
                 'alpha7': ref.select_atoms('resid 168-174 and name BB'),\
                 'alpha8': ref.select_atoms('resid 180-208 and name BB')}

   out='rmsd_'+get_basename(sys.argv[2])[:-4]+'_ref'+get_basename(sys.argv[3])[:-4]+'.txt'
   with open(out, 'w') as f:
      for t in u.trajectory:
         r=rmsd(sel.positions,refsel.positions,center=True,superposition=True)
         f.write('%8.4f %8.4f' % (t.time/1000.0, r))
         for name, s in helices.items():
            r=rmsd(s.positions,refhelices[name].positions,center=True,superposition=True)
            f.write('%8.4f' % r)
         f.write('\n')



if __name__ == '__main__':
   main()

