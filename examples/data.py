#coding:utf8

"""
This is an example of how to use the `PhaseSpace.data` object of p2sat.

It allows to generate and manipulate phase space
"""

# Import p2sat
p2sat_path="../"
import sys
if p2sat_path not in sys.path:sys.path.append(p2sat_path)
import p2sat

# Boolean to export or not the generated phase space
export = False
check_input = False

# Instanciate a PhaseSpace object for electron specie
gps1 = p2sat.PhaseSpace(specie="gamma")

# Define energy and angle parameters
ekin_dict = {"law":"exp","ekin0":1.0}
theta_dict = {"law":"gauss","mu":0.,"sigma":5.}
phi_dict = {"law":"iso"}
x_dict = {"law":"gauss","mu":0.,"sigma":10.}
r_dict = {"law":"gauss","mu":0.,"sigma":50.}
t_dict = {"law":"exp","t0":150.}

# Generate particle phase space
gps1.data.generate(Nconf = 1e4, Npart = 1e12,
                  ekin=ekin_dict, theta=theta_dict, phi=phi_dict,
                  x=x_dict, r=r_dict, t=t_dict)

# Look at the consistency of phase space generation
if check_input:
    print(gps1)
    gps1.plot.figure(0)
    gps1.plot.h1('ekin',bwidth=.1,log=True)
    gps1.plot.f1('ekin',func_name="exp",log=True,bwidth=.1)

    gps1.plot.figure(1)
    gps1.plot.h2('theta','phi',log=True,
                bwidth1=.5,bwidth2=1.)

    gps1.plot.figure(2)
    gps1.plot.h1('r',bwidth=1)
    gps1.plot.f1('r',func_name="gauss",bwidth=1)

# Copy current PhaseSpace in a new object
gps2 = gps1.copy()
# Rotate and translate phase spaces
gps1.data.transformate(T=(-200.,0.,0.),R=(0.,45.,45.))
gps2.data.transformate(T=(200.,0.,0.),R=(0.,0.,180.),rotate_first=True)

# Propagate to a given time
gps1.data.propagate(time=300.)
gps2.data.propagate(time=300.)

# Combine the 2 previous PhaseSpace
gps = gps1 + gps2

# Plot results
gps.plot.figure(3)
gps.plot.h2('x','px',bwidth1=1.,bwidth2=.1,log=True)
gps.plot.figure(4)
gps.plot.h2('x','y',bwidth1=1.,bwidth2=1.,log=True)

# Export phase space if needed
if export:
    gps.export.txt("test_gps.csv",sep=",")
