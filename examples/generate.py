#coding:utf8

"""
This is an example of how to use the `PhaseSpace.stat` object of p2sat.

It allows to make statistics on phase space data in a very simple way
"""

# Import p2sat
p2sat_path="../"
import sys
if p2sat_path not in sys.path:sys.path.append(p2sat_path)
import p2sat

# Boolean to export or not the generated phase space
export = False

# Instanciate a PhaseSpace object for electron specie
eps = p2sat.PhaseSpace(specie="electron")

# Define energy and angle parameters
# ekin_dict = {"law":"mono","E":1.0}
ekin_dict = {"law":"exp","T":1.0}
theta_dict = {"law":"gauss","mu":0.,"sigma":5.}
# theta_dict = {"law":"iso","max":20.}
# theta_dict = {"law":"iso"}
# phi_dict = {"law":"iso","mangle":45.}
phi_dict = {"law":"iso"}

# Generate particle phase space
eps.data.generate(Nconf = 1e5, Npart = 1e12, ekin=ekin_dict, theta=theta_dict, phi=phi_dict)

# Look at the consistency of phase space generation
print(eps)
eps.plot.figure(0)
eps.plot.h1('ekin',bwidth=.1,log=True)
eps.plot.f1('ekin',func_name="exp",log=True,bwidth=.1)

eps.plot.figure(1)
eps.plot.h2('theta','phi',log=True,
            bwidth1=.5,bwidth2=1.)

# Export phase space if needed
if export:
    eps.export.txt("test_eps.csv",sep=",")
