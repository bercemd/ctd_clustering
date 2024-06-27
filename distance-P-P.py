#This python script counts the number of occations where a particular distance between two Phosphorous atoms
#belong to two different peptides and put them into histograms for specific distances. Then average over
#the total number of frames of the trajectory which will give the average count per frame for specific distance
#Modified by Amith, W. D. and Dutagaci, B. 2024, University of California Merced, Merced, CA 95343 USA.
#Based on the examples in the user guide of MDAnalysis package
#Reference:Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations.
#J.Comp.Chem 2011, 32 (10), 2319-2327. DOI: https://doi.org/10.1002/jcc.21787
import MDAnalysis as mda
from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysis.analysis import distances
import numpy as np
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

arg = parser.parse_args()


def dist_mda(u,selecta,selectb):
  sela = u.select_atoms(selecta)
  selb = u.select_atoms(selectb)
 
  sela_seg = sela.segments
  selb_seg = selb.segments

  
  distance = np.zeros(shape=(len(sela), len(selb)))
  dists = []

  for ts in u.trajectory:
    boxsize = ts.dimensions

    for i in range(len(sela_seg)):
        for j in range(len(selb_seg)): 
           if sela_seg[i] != selb_seg[j]:
               rsela = sela.positions
               rselb = selb.positions
               distance = distances.distance_array(rsela,rselb,box=boxsize,backend='OpenMP')
               dists.append(distance)
           else:
               continue
  return dists, sela, selb



     
U = mda.Universe(arg.psf_path, arg.traj_path, in_memory=True)
dists, sela, selb = dist_mda(U,arg.sela,arg.selb) # distance calculation

dist_arr = np.array(dists)


n_frames = 0
nbins = 150

bins = np.zeros(nbins)

boxsize = U.trajectory.ts.dimensions
r_cut = boxsize[0]/2
dr = r_cut/nbins


for ts in U.trajectory:
   
    n_frames = n_frames + 1    

    for i in range(len(sela)): #histogram definition
     for j in range(len(selb)):
         if dist_arr[ts.frame][i][j] != 0.0 and dist_arr[ts.frame][i][j] <= r_cut:  # condition to avoid counting distance between same atoms and go upto half of the box length
             rr = dist_arr[ts.frame][i][j]/dr
             rval = int(np.floor(rr))
             bins[rval] = bins[rval] + 1
         else:
             continue



rvals = np.zeros(nbins)

for i in range(nbins):

    if i==0:
        rvals[i] = 0.5*dr
        bins[i] = bins[i]/n_frames
    else:
        rvals[i] = (i+0.5)*dr
        bins[i] = bins[i]/n_frames

    print(rvals[i], bins[i])


