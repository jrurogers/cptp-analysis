#!/usr/bin/python

"""
calc_rmsd2pdb_allatom.py gro xtc pdb
created by: Julia Rogers
created on: 11-4-19

calculates rmsd of CA using pdb as reference
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
   sel=u.select_atoms('(resid %i-%i) and name CA' % (8+md_resid_offset, 214+md_resid_offset))
   ref=Universe(sys.argv[3])
   refsel=ref.select_atoms('(resid 8-214) and name CA')

   prot = u.select_atoms('protein')
   if prot.resids[0] > 1: md_resid_offset = 264
   else: md_resid_offset = 0


   # helix definitions
   helices = {'alphaN': u.select_atoms('resid %i-%i and name CA' % (10+md_resid_offset, 20+md_resid_offset)),\
              'alpha1': u.select_atoms('resid %i-%i and name CA' % (29+md_resid_offset, 46+md_resid_offset)),\
              'alpha2': u.select_atoms('resid %i-%i and name CA' % (51+md_resid_offset, 69+md_resid_offset)),\
              'alpha3': u.select_atoms('resid %i-%i and name CA' % (79+md_resid_offset, 89+md_resid_offset)),\
              'alpha4': u.select_atoms('resid %i-%i and name CA' % (104+md_resid_offset, 127+md_resid_offset)),\
              'alpha5': u.select_atoms('resid %i-%i and name CA' % (134+md_resid_offset, 144+md_resid_offset)),\
              'alpha6': u.select_atoms('resid %i-%i and name CA' % (152+md_resid_offset, 163+md_resid_offset)),\
              'alpha7': u.select_atoms('resid %i-%i and name CA' % (168+md_resid_offset, 174+md_resid_offset)),\
              'alpha8': u.select_atoms('resid %i-%i and name CA' % (180+md_resid_offset, 208+md_resid_offset))}
   refhelices = {'alphaN': ref.select_atoms('resid 10-20 and name CA'),\
                 'alpha1': ref.select_atoms('resid 29-46 and name CA'),\
                 'alpha2': ref.select_atoms('resid 51-69 and name CA'),\
                 'alpha3': ref.select_atoms('resid 79-89 and name CA'),\
                 'alpha4': ref.select_atoms('resid 104-127 and name CA'),\
                 'alpha5': ref.select_atoms('resid 134-144 and name CA'),\
                 'alpha6': ref.select_atoms('resid 152-163 and name CA'),\
                 'alpha7': ref.select_atoms('resid 168-174 and name CA'),\
                 'alpha8': ref.select_atoms('resid 180-208 and name CA')}

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

