#coding:utf8

"""
This is an example of how to use the `hist` sub-object of PhaseSpace object.

It allows to make histogram from phase space data in a very simple way
"""

# Import p2sat
p2sat_path="../"
import sys
if p2sat_path not in sys.path:sys.path.append(p2sat_path)
import p2sat

# Instanciate a PhaseSpace object for electron specie
eps = p2sat.PhaseSpace(specie="electron")

# Import data from a file
eps.extract.txt("input.tsv",sep=None)

# Get spectrum (Number/MeV, bin width of 0.1 MeV)
ekin,spec = eps.hist.h1('ekin',bwidth=0.1)

# Get spectrum for particles with theta angle between -10 and 10 deg
ekin,spec = eps.hist.h1('ekin', bwidth=0.1, select={'theta':[-10.0,10.0]})

# Fit the last spectrum for ekin > 0.511 MeV (Ne is total number of e-, Te its temperature in MeV)
fekin,Ne,Te = eps.hist.f1('ekin', func_name="exp", bwidth=0.1, select={'theta':[-10.0,10.0],'ekin':[0.511,None]})

print(Ne,Te)
