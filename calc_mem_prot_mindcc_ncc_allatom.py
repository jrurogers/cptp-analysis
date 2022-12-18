#!/usr/bin/python

"""
calc_mem_prot_mindcc_ncc.py gro xtc
created by: Julia Rogers
created on: 9-23-20

calc min dcc, ncc for hydrophobic lipid-prot interactions

"""

NCCCUTOFF = 10.0 # cutoff for close contact, angstrom

import sys, re, os
import numpy as np
from MDAnalysis.analysis.distances import distance_array, self_distance_array
from MDAnalysis.analysis.contacts import soft_cut_q
from MDAnalysis.lib.util import convert_aa_code
from MDAnalysis import *

# hydrophobic carbons of each amino acid
aa_hc = {'ALA': 'CB', 'ILE': 'CB CG2 CG1 CD', 'LEU': 'CB CG CD1 CD2',\
         'MET': 'CB CG SD CE', 'PHE': 'CB CG CD1 CE1 CZ CD2 CE2',\
         'TRP': 'CB CG CD1 CE2 CD2 CE3 CZ3 CZ2 CH2', 'TYR': 'CB CG CD1 CE1 CD2 CE2',\
         'VAL': 'CB CG1 CG2', 'PRO': 'CD CB CG', 'CYS': 'CB'}

def main():

   # select atoms
   u=Universe(sys.argv[1], sys.argv[2])
   prot = u.select_atoms('protein')

   # determine ref leaflet that protein is bound to
   comtop = u.select_atoms('resid 1:132').center_of_mass(pbc=True)[2]
   combot = u.select_atoms('resid 133:264').center_of_mass(pbc=True)[2]
   comprot = prot.center_of_mass(pbc=True)[2]
   dist1 = comprot - comtop
   dist2 = combot - comprot
   box=u.dimensions
   dist1 -= round(dist1/box[2],0)*box[2]
   dist2 -= round(dist2/box[2],0)*box[2]
   if abs(dist1) < abs(dist2): topref = True
   else: topref = False
   if topref: memc = u.select_atoms('resid 1:132 and name C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218 C33 C34 C35 C36 C37 C38 C39 C310 C311 C312 C313 C314 C315 C316 C3F C4F C5F C6F C7F C8F C9F C10F C11F C12F C13F C15F C14F C16F C5S C6S C7S C8S C9S C10S C11S C12S C13S C14S C15S C16S C17S C18S')
   else: memc = u.select_atoms('resid 133:264 and name C23 C24 C25 C26 C27 C28 C29 C210 C211 C212 C213 C214 C215 C216 C217 C218 C33 C34 C35 C36 C37 C38 C39 C310 C311 C312 C313 C314 C315 C316 C3F C4F C5F C6F C7F C8F C9F C10F C11F C12F C13F C15F C14F C16F C5S C6S C7S C8S C9S C10S C11S C12S C13S C14S C15S C16S C17S C18S')

   protsel = ' or '.join(['(resname %s and name %s)' % (aa,aa_hc[aa]) for aa in aa_hc.keys()])
   protc = u.select_atoms(protsel)

   print "Calculating min dcc and ncc between lipids and protein residues for gro %s and xtc %s" % (sys.argv[1], sys.argv[2])
   if topref: print "Using top leaflet as reference"
   else: print "Using bottom leaflet as reference"

   # number of carbons per lipid
   numc_per_lipid = np.empty(len(memc.residues), dtype=int)
   for i in range(len(memc.residues)):
      numc_per_lipid[i] = len(memc.select_atoms('resid %i' % memc.residues[i].resid))
   # number of carbons per protein residue
   numc_per_residue = np.empty(len(protc.residues), dtype=int)
   for i in range(len(protc.residues)):
      numc_per_residue[i] = len(protc.select_atoms('resid %i' % protc.residues[i].resid))

   # loop over traj and calc quantities
   fnamel = 'memlip_prot_mindcc_ncc_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   fnamep = 'prot_mindcc_ncc_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   with open(fnamel, 'w') as fl, open(fnamep, 'w') as fp:
      for t in u.trajectory:
         if t.time/1000.0 < 50: continue # statistics from 50-100 ns
         if t.time/1000.0 % 10 == 0: print t.time/1000.0
 
         # distances between all hydrophobic carbons of lipids and protein
         dccs = distance_array(memc.positions, protc.positions, box=t.dimensions)

         # for each lipid
         for i in range(len(memc.residues)):
            # get all relevant carbon-carbon dists
            sndx = np.sum(numc_per_lipid[:i])
            endx = np.sum(numc_per_lipid[:i+1])
            # calc min dcc and ncc
            mindcc = np.amin(dccs[sndx:endx,:])
            ncc = len(np.where(dccs[sndx:endx,:] <= NCCCUTOFF)[0])
            fl.write('%8.4f %5i %8.4f %8i\n' % (t.time/1000.0,memc.residues[i].resid,mindcc,ncc))

         # for each protein residue
         for i in range(len(protc.residues)):
            # get all relevant carbon-carbon dists
            sndx = np.sum(numc_per_residue[:i])
            endx = np.sum(numc_per_residue[:i+1])
            # calc min dcc and ncc
            mindcc = np.amin(dccs[:,sndx:endx])
            ncc = len(np.where(dccs[:,sndx:endx] <= NCCCUTOFF)[0])
            fp.write('%8.4f %5i %8.4f %8i\n' % (t.time/1000.0,protc.residues[i].resid-264,mindcc,ncc))



if __name__ == '__main__':
   main()

