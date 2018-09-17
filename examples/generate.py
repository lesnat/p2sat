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
ekin_dict = {"law":"exp","ekin0":1.0}
theta_dict = {"law":"gauss","mu":0.,"sigma":5.}
phi_dict = {"law":"iso"}
x_dict = {"law":"mono","x0":1.5}
r_dict = {"law":"gauss","mu":0.,"sigma":50.}
t_dict = {"law":"exp","t0":150.}

# Generate particle phase space
eps.data.generate(Nconf = 1e4, Npart = 1e12,
                  ekin=ekin_dict, theta=theta_dict, phi=phi_dict,
                  x=x_dict, r=r_dict, t=t_dict)

# Look at the consistency of phase space generation
print(eps)
eps.plot.figure(0)
eps.plot.h1('ekin',bwidth=.1,log=True)
eps.plot.f1('ekin',func_name="exp",log=True,bwidth=.1)

eps.plot.figure(1)
eps.plot.h2('theta','phi',log=True,
            bwidth1=.5,bwidth2=1.)

eps.plot.figure(2)
eps.plot.h1('r',bwidth=.1)
eps.plot.f1('r',func_name="gauss",bwidth=.1)

# Export phase space if needed
if export:
    eps.export.txt("test_eps.csv",sep=",")
