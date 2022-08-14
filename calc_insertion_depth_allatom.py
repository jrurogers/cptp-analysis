#!/usr/bin/python

"""
calc_insertion_depth_allatom.py gro xtc
created by: Julia Rogers
created on: 6-24-20

calculate average distance in z between each protein residue COM
and average phosphate position

** for membrane bound configurations, uses first frame to determine ref leaflet **

"""


import sys
import numpy as np
from MDAnalysis import *
from math import floor


def com_pbc(mol, box):
   com = np.zeros(3)
   for i in range(len(mol)):
      for j in range(3):
         com[j] += mol[i][j] + int((mol[0][j]-mol[i][j])/(box[j]*0.5))*box[j]
   com /= len(mol)
   return com

def main():
   u=Universe(sys.argv[1], sys.argv[2])
   u.trajectory[0]
   top=u.select_atoms('resid 1:132')
   bot=u.select_atoms('resid 133:264')
   prot=u.select_atoms('protein and not name H*')
   box=u.dimensions

   print(top.residues)
   print(bot.residues)
   print(prot.residues)

   # determine ref leaflet
   comtop = top.center_of_mass(pbc=True)[2]
   combot = bot.center_of_mass(pbc=True)[2]
   comprot = prot.center_of_mass(pbc=True)[2]
   dist1 = comprot - comtop
   dist2 = combot - comprot
   dist1 -= round(dist1/box[2],0)*box[2]
   dist2 -= round(dist2/box[2],0)*box[2]
   if abs(dist1) < abs(dist2): topref = True
   else: topref = False
   print(dist1, dist2, topref)

   # phosphate selection
   if topref:
      po4 = top.select_atoms('name P')
   else:
      po4 = bot.select_atoms('name P')

   numres = len(prot.residues)
   depths = [[] for i in range(numres)]
   for t in u.trajectory:
      if t.time/1000.0 % 10 == 0: print(t.time/1000.0)
      zref = np.mean(po4.positions[:,2])
      for i in range(numres):
         res = prot.select_atoms('resid %i' % (i+265))
         if topref: depth = res.center_of_mass(pbc=True)[2] - zref
         else: depth = zref - res.center_of_mass(pbc=True)[2]
         depths[i].append(depth)

   fname='av_prot_insertion_depth.txt'
   with open(fname, 'w') as f:
      for i in range(numres):
         f.write('%5i %8.4f\n' % ((i+1), np.mean(depths[i])))


if __name__ == '__main__':
   main()

