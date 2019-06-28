#coding:utf8

"""
This is an example of how to use the `PhaseSpace.data` object of p2sat.

It allows to generate and manipulate phase space
"""

# Import p2sat
p2sat_path="../"
import sys
if p2sat_path not in sys.path:sys.path.insert(0,p2sat_path)
import p2sat

# Boolean to export or not the generated phase space
export = False
check_input = True

# Instanciate a PhaseSpace object for electron specie
gps1 = p2sat.PhaseSpace(particle="gamma")

# Define energy and angle parameters
ekin_dict   = dict(law="exp",scale=1.)                  # exponential on energy
phi_dict    = dict(law="uni",min=0,max=360)             # isotropic on phi
omega_dict  = dict(law="gauss",mu=0,sigma=5.)           # gaussian on solid angle
x_dict      = dict(law="gauss",mu=0,sigma=10.)          # gaussian on x
r_dict      = dict(law="gauss",mu=0,sigma=10.)          # gaussian on r
t_dict      = dict(law="exp",scale=150.)                # exponential on time

# Generate particle phase space
gps1.data.generate(Nconf = 1e4, Npart = 1e12,
                  ekin=ekin_dict, phi=phi_dict, omega=omega_dict,
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
gps1.data.transformate(T=(-200.,0.,0.),R=(0.,0.,45.))
gps2.data.transformate(T=(200.,0.,0.),R=(0.,0.,180.),rotate_first=True)

# Propagate to a given time
gps1.data.propagate(t=300.)
gps2.data.propagate(t=300.)

# Combine the 2 previous PhaseSpace
gps = gps1 + gps2

# Plot results
if check_input:
    gps.plot.figure(3)
    gps.plot.h2('x','px',bwidth1=1.,bwidth2=.1,log=True)
    gps.plot.figure(4)
    gps.plot.h2('x','y',bwidth1=1.,bwidth2=1.,log=True)
    gps.plot.figure(5)
    gps.plot.h3('x','y','z',wmin=1e6,s=10,log=True)

# Export phase space if needed
if export:
    # Discretize phase space, to limitate disk usage
    bwidth=[1.,2.,2.,1.,1.,1.]
    bwidth = None
    gps.data.discretize(with_time=False,split=2,MP=False,bwidth=bwidth)
    gps.export.txt("test_gps.csv",sep=",")
