#This python script calculates the self diffusion coefficient from Einstein formula for a specific
#selection by linear regression of mean-squared displacment (MSD) values at given lag-times within a specified time range
#Modified by Amith, W. D. and Dutagaci, B. 2024, University of California Merced, Merced, CA 95343 USA.
#Based on the examples in the user guide of MDAnalysis package
#Reference:Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. MDAnalysis: A toolkit for the analysis of molecular dynamics simulations.
#J.Comp.Chem 2011, 32 (10), 2319-2327. DOI: https://doi.org/10.1002/jcc.21787
import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
from MDAnalysis.tests.datafiles import PSF, DCD
import numpy as np
from scipy.stats import linregress


u = mda.Universe("step4_input.psf", "unwrapped_production.dcd")
MSD = msd.EinsteinMSD(u, select='protein', msd_type='xyz', fft=True)
MSD.run()

msd =  MSD.results.timeseries

nframes = MSD.n_frames
timestep = 10 # this needs to be the actual time between frames
lagtimes = np.arange(nframes)*timestep # make the lag-time axis

start_time = 20
start_index = int(start_time/timestep)
end_time = 100
end_index = int(end_time/timestep)
linear_model = linregress(lagtimes[start_index:end_index],
                                              msd[start_index:end_index])
slope = linear_model.slope
# dim_fac is 3 as we computed a 3D msd with 'xyz'
D = slope * 1/(2*MSD.dim_fac)
print ("D =", D)
