#This python script calculates the average number of normalized contacts below or equal to 5 Å for each residue pair 
#between different peptides over the total number of combinations for each residue pair
#Modified by Amith, W. D. and Dutagaci, B. 2024, University of California Merced, Merced, CA 95343 USA.
#Based on the examples in the user guide of MDAnalysis package
#Reference:Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations.
#J.Comp.Chem 2011, 32 (10), 2319-2327. DOI: https://doi.org/10.1002/jcc.21787

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.contacts import contact_matrix
from MDAnalysis.analysis.distances import distance_array
from pylab import *
import matplotlib.pyplot as plt

import os,sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '--traj_path',type=str, default='md.dcd',
    help='Trajectory path.')
parser.add_argument(
    '--psf_path',type=str, default='protein.psf',
    help='PSF path.')
parser.add_argument(
    '--sela',type=str, default='protein',
    help='First selection.')
parser.add_argument(
    '--selb',type=str, default='segid DNAA',
    help='Second selection.')
parser.add_argument(
    '--out_path',type=str, default='rmsf.dat',
    help='Output path.')
arg = parser.parse_args()

def contacts_mda(u,selecta,selectb):
        #initial selection of desired atoms both whole subunit and atoms close to target comparison
        sel = u.select_atoms(selecta, updating=True)
        ref = u.select_atoms(selectb, updating=True)

        #getting residue names and residue number from selected groups
        sel_res = sel.residues
        sel_resids = sel_res.resids
        sel_nums = sel_res.resnums
        sel_segnames = sel_res.segids
        ref_res = ref.residues
        ref_resids = ref_res.resids
        ref_nums = ref_res.resnums
        ref_segnames = ref_res.segids
        #empty array to stack all frame's counts of contacts at positions i,j
        count_stack = np.zeros((len(sel_res), len(ref_res)))

        #iteration over the trajectory
        for frame in u.trajectory:

                #array of 100's to be replaced at positions i,j by the minimum distance between residue i and residue j
                res_dist = np.full((len(sel_res), len(ref_res)), 100)

                #iterating through each of the contacting residues
                for n in sel_res:
                        for m in ref_res:
                                #selecting the atoms of sel and ref residues
                                sel_atoms = n.atoms.positions
                                ref_atoms = m.atoms.positions

                                #calculating the minimum distance between the residues
                                distance = distance_array(sel_atoms, ref_atoms, backend='OpenMP')
                                min_dist = np.min(distance)

                                #replacing the minimum distance value of residue i and residue j at position [i,j] in res_dist
                                sel_i = (np.where(n == sel_res)[0][0])
                                ref_j = (np.where(m == ref_res)[0][0])
                                res_dist[sel_i, ref_j] = min_dist
                #Calculating where in res_dist is less than or equal to a cutoff distance
                cutoff=5
                contacts = contact_matrix(res_dist, cutoff)
                #Converting True to 1, and False to 0
                contact_int = contacts * 1
                #adding the contacts in the frame to the total contacts in count_stack
                count_stack = count_stack + contact_int

        #Averaging the total number of contacts by the number of frames in the trajectory (needs to be transposed after division)
        count_avg = count_stack / len(u.trajectory)
        #counts = count_avg.transpose()
        counts = count_avg

        return counts, sel_resids, sel_segnames, ref_resids, ref_segnames


U = mda.Universe(arg.psf_path, arg.traj_path, in_memory=True)
contacts_arr, sela_resids, sela_segnames, selb_resids, selb_segnames = contacts_mda(U,arg.sela,arg.selb)
contacts_arr.shape
outputfile = open(arg.out_path,"w")
for i in range(len(contacts_arr)):
  for j in range(len(contacts_arr[i])):
      if sela_segnames[i] != selb_segnames[j] and i < j:
          print(sela_resids[i],selb_resids[j],contacts_arr[i][j],sela_segnames[i],selb_segnames[j],file=outputfile)
      else:
          continue
outputfile.close()


resa = []
resb = []
contacts = []

inputfile = open(arg.out_path,"r")
value = inputfile.readlines()
inputfile.close()

for line in value:
    resa.append(int(line.split()[0]))
    resb.append(int(line.split()[1]))
    contacts.append(float(line.split()[2]))

resa = np.asarray(resa)
resb = np.asarray(resb)
contacts = np.asarray(contacts)

for i in range(1,15):
    for j in range(1,15):
        ctemp = []
        for k in range(len(resa)):
         if i == resa[k] and j == resb[k]:
             ctemp.append(contacts[k])
        print(i, j, len(ctemp), sum(ctemp)/len(ctemp))
