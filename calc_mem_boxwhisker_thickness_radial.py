#!/usr/bin/python

"""
calc_mem_boxwhisker_thickness_radial.py gro xtc
created by: Julia Rogers
created on: 11-3-22

calc box whisker of thickness based on radial distance
(temporal avs)


*** MUST CENTER PROT IN BOX AND CORRECT FOR PBC FIRST ***
*** e.g. gmx trjconv -pbc res -ur tric -ceneter       ***
*** cannot properly account for pbc once rotated      ***

"""

REFPROT = "reference_structures/av_xy_plane_ref.gro"

import sys, os
import numpy as np
import math
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
   if topref:
      # top leaflet taken to be one that CPTP is bound to
      topmemsel = 'resid 1:132 and name C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218 C33 C34 C35 C36 C37 C38 C39 C310 C311 C312 C313 C314 C315 C316 C3F C4F C5F C6F C7F C8F C9F C10F C11F C12F C13F C15F C14F C16F C5S C6S C7S C8S C9S C10S C11S C12S C13S C14S C15S C16S C17S C18S'
      botmemsel = 'resid 133:264 and name C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218 C33 C34 C35 C36 C37 C38 C39 C310 C311 C312 C313 C314 C315 C316 C3F C4F C5F C6F C7F C8F C9F C10F C11F C12F C13F C15F C14F C16F C5S C6S C7S C8S C9S C10S C11S C12S C13S C14S C15S C16S C17S C18S'
   else:
      topmemsel = 'resid 133:264 and name C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218 C33 C34 C35 C36 C37 C38 C39 C310 C311 C312 C313 C314 C315 C316 C3F C4F C5F C6F C7F C8F C9F C10F C11F C12F C13F C15F C14F C16F C5S C6S C7S C8S C9S C10S C11S C12S C13S C14S C15S C16S C17S C18S'
      botmemsel = 'resid 1:132 and name C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218 C33 C34 C35 C36 C37 C38 C39 C310 C311 C312 C313 C314 C315 C316 C3F C4F C5F C6F C7F C8F C9F C10F C11F C12F C13F C15F C14F C16F C5S C6S C7S C8S C9S C10S C11S C12S C13S C14S C15S C16S C17S C18S'

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

   print "Creating box whisker of thickness for gro %s and xtc %s" % (sys.argv[1], sys.argv[2])
   if topref: print "Using top leaflet as reference"
   else: print "Using bottom leaflet as reference"

   # distance grid
   dgrid = np.linspace(0.0,45.0,num=10)
   dbin_thickness = [[] for i in range(len(dgrid)-1)]

   topmemc = u.select_atoms(topmemsel)
   botmemc = u.select_atoms(botmemsel)
   for t in u.trajectory:
      if t.time/1000.0 < 50: continue # statistics from 50-100 ns
      if t.time/1000.0 % 10 == 0: print t.time/1000.0

      # positions of com of lipids
      topxs = np.zeros(len(topmemc.residues))
      topys = np.zeros(len(topmemc.residues))
      botxs = np.zeros(len(botmemc.residues))
      botys = np.zeros(len(botmemc.residues))
      # z coordinate of lipid's phosphate
      topphos = np.zeros(len(topmemc.residues))
      botphos = np.zeros(len(botmemc.residues))
      
      # com of each lipid in top leaflet
      for i in range(len(topmemc.residues)):
         comlip = u.select_atoms('resid %i' % topmemc.residues[i].resid).center_of_mass()[:2]
         p = u.select_atoms('resid %i and name P' % topmemc.residues[i].resid).positions[0,2]
         topxs[i] = comlip[0]
         topys[i] = comlip[1]
         topphos[i] = p

      # com of each lipid in bot leaflet
      for i in range(len(botmemc.residues)):
         comlip = u.select_atoms('resid %i' % botmemc.residues[i].resid).center_of_mass()[:2]
         p = u.select_atoms('resid %i and name P' % botmemc.residues[i].resid).positions[0,2]
         botxs[i] = comlip[0]
         botys[i] = comlip[1]
         botphos[i] = p

      # radial distance of lipids from center
      toprs = np.sqrt(np.power(topxs - boxcenter[0], 2) + np.power(topys - boxcenter[1], 2))
      botrs = np.sqrt(np.power(botxs - boxcenter[0], 2) + np.power(botys - boxcenter[1], 2))
      for i in range(1, len(dgrid)):
         topndxs = (toprs >= dgrid[i-1]) & (toprs < dgrid[i])
         botndxs = (botrs >= dgrid[i-1]) & (botrs < dgrid[i])
         thickness = np.mean(topphos[topndxs]) - np.mean(botphos[botndxs])
         #print(t.time/1000.0, i, dgrid[i], sum(topndxs), sum(botndxs), thickness)
         if not math.isnan(thickness): dbin_thickness[i-1].append(thickness)

   fname = 'memlip_boxwhisker_thickness_radial_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for i in range(1, len(dgrid)):
         f.write('%8.4f' % dgrid[i])
         stats = calc_stats(dbin_thickness[i-1])
         print(len(dbin_thickness[i-1]), stats)
         for j in range(len(stats)):
            f.write('%14.6f' % stats[j])
         f.write('\n')


if __name__ == '__main__':
   main()

