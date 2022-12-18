#!/usr/bin/python

"""
calc_mem_av_thickness_xy.py gro xtc
created by: Julia Rogers
created on: 7-16-20

calc av membrane thickness (z distance btwn phos) in xy plane around protein


*** MUST CENTER PROT IN BOX AND CORRECT FOR PBC FIRST ***
*** e.g. gmx trjconv -pbc res -ur tric -ceneter       ***
*** cannot properly account for pbc once rotated      ***

"""

REFPROT = "reference_structures/av_xy_plane_ref.gro"

import sys, re, os
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

   print "Calculating av thickness on xy grid for gro %s and xtc %s" % (sys.argv[1], sys.argv[2])
   if topref: print "Using top leaflet as reference"
   else: print "Using bottom leaflet as reference"

   # grid
   xgrid = np.linspace(1.0,93.0,num=47)
   ygrid = np.linspace(1.0,93.0,num=47)
   grid_width = 2.0

   # positions of com of lipids
   topxs = []
   topys = []
   botxs = []
   botys = []
   # z coordinate of lipid's phosphate
   topphos = []
   botphos = []

   topmemc = u.select_atoms(topmemsel)
   botmemc = u.select_atoms(botmemsel)
   for t in u.trajectory:
      if t.time/1000.0 < 50: continue
      if t.time/1000.0 % 10 == 0: print t.time/1000.0

      # com of each lipid in top leaflet
      for i in range(len(topmemc.residues)):
         comlip = u.select_atoms('resid %i' % topmemc.residues[i].resid).center_of_mass()[:2]
         p = u.select_atoms('resid %i and name P' % topmemc.residues[i].resid).positions[0,2]
         topxs.append(comlip[0])
         topys.append(comlip[1])
         topphos.append(p)

      # com of each lipid in bot leaflet
      for i in range(len(botmemc.residues)):
         comlip = u.select_atoms('resid %i' % botmemc.residues[i].resid).center_of_mass()[:2]
         p = u.select_atoms('resid %i and name P' % botmemc.residues[i].resid).positions[0,2]
         botxs.append(comlip[0])
         botys.append(comlip[1])
         botphos.append(p)

   topxs = np.array(topxs)
   topys = np.array(topys)
   topphos = np.array(topphos)
   botxs = np.array(botxs)
   botys = np.array(botys)
   botphos = np.array(botphos)
   print topxs.shape, topys.shape, topphos.shape, botxs.shape, botys.shape, botphos.shape

   fname = 'memlip_av_thickness_xy_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fname, 'w') as f:
      for i in range(len(xgrid)):
         top_xndxs = (topxs >= xgrid[i]-0.5*grid_width) & (topxs < xgrid[i]+0.5*grid_width)
         bot_xndxs = (botxs >= xgrid[i]-0.5*grid_width) & (botxs < xgrid[i]+0.5*grid_width)
         for j in range(len(ygrid)):
            top_yndxs = (topys >= ygrid[j]-0.5*grid_width) & (topys < ygrid[j]+0.5*grid_width)
            bot_yndxs = (botys >= ygrid[j]-0.5*grid_width) & (botys < ygrid[j]+0.5*grid_width)
            
            avtopphos = np.mean(topphos[top_xndxs & top_yndxs])
            avbotphos = np.mean(botphos[bot_xndxs & bot_yndxs])
            thickness = avtopphos - avbotphos
            f.write('%8.4f %8.4f %8.4f\n' % (xgrid[i],ygrid[j],thickness))
         f.write('\n')
   


if __name__ == '__main__':
   main()

