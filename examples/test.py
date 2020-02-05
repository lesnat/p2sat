#coding:utf8

"""
...
"""

# Import p2sat
p2sat_path="../"
import sys
if p2sat_path not in sys.path:sys.path.insert(0,p2sat_path)
import p2sat

eps = p2sat.datasets.PhaseSpace(specie="e-", unit_system="UHI")
eps.load.txt("ExamplePhaseSpace.csv", in_code_units=False, sep=",")

p2sat.plot.figure(0, clear=True)
p2sat.plot.hist1d(eps, "ekin", legend="test $E_{kin}$", log=True)
