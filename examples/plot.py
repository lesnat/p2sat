#coding:utf8

"""
This is an example of how to use the `PhaseSpace.plot` object of p2sat.

It allows to make plots from phase space data in a very simple way
"""

# Import p2sat
p2sat_path="../"
import sys
if p2sat_path not in sys.path:sys.path.append(p2sat_path)
import p2sat

# Instanciate a PhaseSpace object for electron specie
eps = p2sat.PhaseSpace(particle="electron")

# Import data from a file
eps.extract.txt("example.csv",sep=",",verbose=False)

# Plot spectrum in log scale (Number/MeV, bin width of 0.1 MeV)
eps.plot.figure(0)
eps.plot.h1('ekin', log=True, bwidth=0.1)

# Fit the last spectrum for ekin > 0.511 MeV
eps.plot.f1('ekin', func_name="exp", log=True, bwidth=0.1, select={'ekin':[0.511,None]})

# Plot Transverse particle dispersion for electrons with kinetic energy > 0.511 MeV (bin width of 10 µm between -500 and 500 µm)
eps.plot.figure(1)
eps.plot.h2('y','z',log=True,
            bwidth1=5.0,bwidth2=5.0,
            brange1=[-300.,300.],brange2=[-300.,300.],
            select={'x':300,'ekin':[0.511,None]})

# Add a contour plot
eps.plot.c2('y','z',log=True, gfilter=3.5,
            bwidth1=5.0,bwidth2=5.0,
            brange1=[-300.,300.],brange2=[-300.,300.],
            select={'x':300,'ekin':[0.511,None]})

# Plot angle/energy polar distribution of the particles
eps.plot.figure(2)
eps.plot.h2('theta','ekin',
            log=True,polar=True,
            bwidth1=1.0,bwidth2=0.1)
