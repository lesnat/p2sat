#coding:utf8
"""
Problematic
===========

When you deal with phase space analysis, the philosophy is often the same :

- Import data from a file or generate it with defined laws
- Make calculations with this raw data (kinetic energy, divergence, ...)
- Filter some of it (divergence with a condition on energy, ...)
- Make histograms
- Make fits
- Plot results
- Format the axes

It is often long and leads to many programmation errors, and when there is a need
to compare several set of datas, it become more and more complicated.

p2sat objective is to unify a set of phase space data in a single object, and
give all the needed methods to manipulate and analyse these datas.

This way the analysis and comparison of several sources need much less efforts
and physicists can do what they like to do : physics.

Package Structure
=================

All the package structure relies on one single object : `PhaseSpace`.

This object is a "box" containing all the informations about a given particle
phase space, and all the methods to interact with it (make histograms, statistics or plots).

Once you create such object, the data it contains are isolated from the
rest of your script (see examples below to show how to acess data), so there is no
risk of confuse between different data sources.
This object-oriented approach is so completely coherent when you need to compare
data from several simulation files for examples : the methods you use are the same,
only the instances names are changing.

Particle phase space
====================

The phase space of a set of particles gives the most complete description of a
dynamical system (assuming particles are independant), and all the needed informations
can be reconstructed from it.

This phase space contains informations such as particles positions, momentum and time.
There is also given a statistical weight to each configuration.
The phase space is then a 7D space (3 for momentum, 3 for space, 1 for time).
If it is discretized, this represent a tremendous amount of possible configurations,
so give a statistical weight to **EACH** of them is not the good approach.
The p2sat approach is to save only the informations on configurations that have
a non-zero statistical weight, and list them, one configuration per line.
To have access to a given configuration, there is only the need of knowing its index.

Informations about units can be found in the `_Data` object (access via `PhaseSpace.data`)

Quick example
=============

Assuming `p2sat` is imported, you can instanciate a `PhaseSpace` object for,
let say, electrons, and import a simulation file containing the phase space
informations.

>>> eps = p2sat.PhaseSpace(particle="electron")
>>> eps.extract.txt("example.csv",sep=",")

All the data in you simulation file can now be found at `eps.data`

>>> # List of all the statistical weights
>>> print(eps.data.raw.w)
[1456.0, 1233.0 , 756.0, ... ]
>>> # List of all the x position
>>> print(eps.data.raw.x)
[10.0, 50.0, 30.0, ... ]
>>> # List of all the kinetic energies
>>> print(eps.data.raw.ekin)
[1.58, 4.61, 3.28, ... ]

This means that at the first index, there is an electron at position
x=10.0 um with 1.58 MeV of kinetic energy and with statistical weight of 1456.0.

You can now do statistics with this data, for example get the standard deviation
of theta angle for all the electrons

>>> theta_std = eps.stat.standard_deviation('theta')

and get an histogram (Number/degree, bin width = 1 degree) of this quantity

>>> theta,Ntheta = eps.hist.h1('theta', bwidth=0.1)

It is also possible to make simple or complicated plot in a elegant way

>>> eps.plot.figure(0)
>>> eps.plot.h1('theta', log=True, bwidth=0.1)
>>> eps.plot.figure(1)
>>> eps.plot.h2('theta','ekin',
...             log=True, polar=True
...             bwidth1=1.0,bwidth2=0.1,
...             select={'x':10.0,'t':[0.0,100.0]})

Other examples and a more complete documentation can be found at :
https://github.com/lesnat/p2sat

"""
__version__ = "1.4.2"

from .PhaseSpace import PhaseSpace

def testing(name="all"):
    from ._tests import run
    if name=="all":
        run(name="raw")
        run(name="hist")
        run(name="plot")
        run(name="heritage")
    else:
        run(name=name)
