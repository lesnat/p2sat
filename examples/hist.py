#coding:utf8

"""
This is an example of how to use the `PhaseSpace.hist` object of p2sat.

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

# Fit the last spectrum for ekin > 0.511 MeV (Ne is total number of e-, Te its temperature in MeV)
fekin,Ne,Te = eps.hist.f1('ekin', func_name="exp", bwidth=0.1, select={'ekin':[0.511,None]})
print("Hot electron temperature (fit) : %.3E MeV"%Te)
