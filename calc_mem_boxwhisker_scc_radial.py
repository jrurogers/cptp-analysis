#!/usr/bin/python

"""
calc_mem_boxwhisker_scc_radial.py gro xtc memlip_scc.txt
created by: Julia Rogers
created on: 11-3-22

calc box whisker of coarse-grained order param Scc based on radial distance

* note: must have already calculated Scc for each lipid *

*** MUST CENTER PROT IN BOX AND CORRECT FOR PBC FIRST ***
*** e.g. gmx trjconv -pbc res -ur tric -ceneter       ***
*** cannot properly account for pbc once rotated      ***

"""

REFPROT = "reference_structures/av_xy_plane_ref.gro"

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

def calc_stats(values):
   q = np.quantile(values, [0, 0.25, 0.5, 0.75, 1])
   return np.mean(values), np.std(values), q[0], q[1], q[2], q[3], q[4] 

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
   if topref: memsel = 'resid 1:132 and name C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218 C33 C34 C35 C36 C37 C38 C39 C310 C311 C312 C313 C314 C315 C316 C3F C4F C5F C6F C7F C8F C9F C10F C11F C12F C13F C15F C14F C16F C5S C6S C7S C8S C9S C10S C11S C12S C13S C14S C15S C16S C17S C18S'
   else: memsel = 'resid 133:264 and name C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218 C33 C34 C35 C36 C37 C38 C39 C310 C311 C312 C313 C314 C315 C316 C3F C4F C5F C6F C7F C8F C9F C10F C11F C12F C13F C15F C14F C16F C5S C6S C7S C8S C9S C10S C11S C12S C13S C14S C15S C16S C17S C18S'

   # ref prot gro
   ref = Universe(REFPROT)
   # shift to box center
   com = ref.select_atoms('protein').center_of_mass(pbc=True)
   boxcenter = ref.dimensions*0.5
   print(boxcenter)
   shift = boxcenter[:3] - com
   ref.trajectory.ts.positions += np.array([shift[0], shift[1], 0])
   ag = ref.atoms
   # orient so along x-axis
   protvec = ref.select_atoms('resid 352 and name CA').positions[0] - ref.select_atoms('resid 445 and name CA').positions[0]
   angle = np.degrees(np.arccos(protvec[0] / np.linalg.norm(protvec)))
   transformations.rotate.rotateby(-angle, direction=[0,0,1], ag=ag)(ref.trajectory.ts)
   with Writer("rotated.gro", ref.atoms.n_atoms) as W:
      W.write(ref.atoms)

   # align traj to ref gro
   workflow=[]
   if not topref:
      ag=u.atoms
      workflow.append(transformations.rotate.rotateby(180.0, direction=[1,0,0], ag=ag))
   # remove trans+rot in xy
   workflow.append(transformations.fit_rot_trans(prot, ref, plane="xy", weights="mass"))
   u.trajectory.add_transformations(*workflow)

   # get scc data
   sccs = np.loadtxt(sys.argv[3], usecols=(2,))

   print "Creating box whisker of APL for gro %s and xtc %s" % (sys.argv[1], sys.argv[2])
   if topref: print "Using top leaflet as reference"
   else: print "Using bottom leaflet as reference"

   # positions of com of lipids
   xs = []
   ys = []
   memc = u.select_atoms(memsel)
   for t in u.trajectory:
      if t.time/1000.0 < 50: continue # statistics from 50-100 ns
      if t.time/1000.0 % 10 == 0: print t.time/1000.0

      # com of each lipid
      for i in range(len(memc.residues)):
         comlip = u.select_atoms('resid %i' % memc.residues[i].resid).center_of_mass()[:2]
         xs.append(comlip[0])
         ys.append(comlip[1])

   xs = np.array(xs)
   ys = np.array(ys)

   # radial distance of lipids from center
   rs = np.sqrt(np.power(xs - boxcenter[0], 2) + np.power(ys - boxcenter[1], 2))
   # distance grid
   dgrid = np.linspace(0.0,45.0,num=10)


   fname = 'memlip_boxwhisker_scc_radial_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for i in range(1, len(dgrid)):
         f.write('%8.4f' % dgrid[i])
         ndxs = (rs >= dgrid[i-1]) & (rs < dgrid[i])
         stats = calc_stats(sccs[ndxs].flatten())
         for j in range(len(stats)):
            f.write('%14.6f' % stats[j])
         f.write('\n')


if __name__ == '__main__':
   main()

