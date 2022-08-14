#!/usr/bin/python

"""
calc_rmsf_martini.py gro xtc
created by: Julia Rogers
created on: 3-31-20

calculates rmsf of backbone beads

"""


import sys
import numpy as np
from MDAnalysis import *
from MDAnalysis.analysis.rms import RMSF
from helper.general import get_basename


def main():
   u=Universe(sys.argv[1], sys.argv[2])
   sel=u.select_atoms('name BB')

   rmsfer = RMSF(sel).run()
   fname='rmsf_'+get_basename(sys.argv[2])[:-4]+'.txt'   
   with open(fname, 'w') as f:
      for i in range(len(rmsfer.rmsf)):
         f.write('%5i %8.4f\n' % ((i+1), rmsfer.rmsf[i]))


if __name__ == '__main__':
   main()

