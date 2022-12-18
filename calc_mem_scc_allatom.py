#!/usr/bin/python

"""
calc_mem_scc_allatom.py gro xtc
created by: Julia Rogers
created on: 11-3-22

calculate coarse-grained order parameter
see: https://lipyphilic.readthedocs.io/en/latest/reference/lib/order_parameter.html
and: https://pubs.acs.org/doi/full/10.1021/acs.jpclett.0c01317

"""


import sys, os
import numpy as np
from MDAnalysis import *
from lipyphilic.lib.order_parameter import SCC


def main():

   # select atoms
   u=Universe(sys.argv[1], sys.argv[2])
   prot = u.select_atoms('protein')
   times = [t.time for t in u.trajectory]

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

   if topref: resid_sel='1:132'
   else: resid_sel='133:264'

   print("Calculating coarse-grained order param for gro %s and xtc %s" % (sys.argv[1], sys.argv[2]))
   if topref: print("Using top leaflet as reference")
   else: print("Using bottom leaflet as reference")

   names = u.select_atoms("name C2? C2?? C?F C??F and not name C2S").names
   for n in set(names):
      print(n)
   names = u.select_atoms("name C3? C3?? C?S C??S and not name C3F C1S C2S").names
   for n in set(names):
      print(n)

   # sn1/acyl tails order param
   tail_sel = "resid %s and name C2? C2?? C?F C??F and not name C2S" % resid_sel
   scc_sn1 = SCC(
      universe=u,
      tail_sel=tail_sel
   )
   scc_sn1.run(start=None, stop=None, step=None, verbose=True)

   # sn2/sphingoid tail order param
   tail_sel = "resid %s and name C3? C3?? C?S C??S and not name C3F C1S C2S" % resid_sel
   scc_sn2 = SCC(
      universe=u,
      tail_sel=tail_sel
   )
   scc_sn2.run(start=None, stop=None, step=None, verbose=True)

   # weighted av
   scc = SCC.weighted_average(scc_sn1, scc_sn2)

   fname = 'memlip_scc_'+os.path.basename(sys.argv[2])[:-4]+'.txt'
   memresidues = u.select_atoms('resid %s' % resid_sel).residues
   nlip = scc.SCC.shape[0]
   with open(fname, 'w') as f:
      for j, t in enumerate(times):
         for i in range(nlip):
            f.write('%8.4f %5i %8.4f %8.4f %8.4f\n' % (t/1000.0, memresidues[i].resid, scc.SCC[i,j], scc_sn1.SCC[i,j], scc_sn2.SCC[i,j]))



if __name__ == '__main__':
   main()

