#!/usr/bin/python

"""
calc_rmsf_allatom.py gro xtc
created by: Julia Rogers
created on: 3-31-20

calculates rmsf of Calphas

"""


import sys, os
import numpy as np
from MDAnalysis import *
from MDAnalysis.analysis.rms import RMSF


def main():
   u=Universe(sys.argv[1], sys.argv[2])
   sel=u.select_atoms('name CA')

   rmsfer = RMSF(sel).run()
   fname='rmsf_'+os.path.basename(sys.argv[2])[:-4]+'.txt'   
   with open(fname, 'w') as f:
      for i in range(len(rmsfer.rmsf)):
         f.write('%5i %8.4f\n' % ((i+1), rmsfer.rmsf[i]))


if __name__ == '__main__':
   main()

