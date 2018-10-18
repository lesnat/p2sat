#coding:utf8

"""
This is an example of how to use the `PhaseSpace.stat` object of p2sat.

It allows to make statistics on phase space data in a very simple way
"""

# Import p2sat
p2sat_path="../"
import sys
if p2sat_path not in sys.path:sys.path.insert(0,p2sat_path)
import p2sat

# Instanciate a PhaseSpace object for electron specie
eps = p2sat.PhaseSpace(particle="electron")

# Import data from a file
eps.extract.txt("example.csv",sep=",",verbose=False)

# Get the mean value of theta
theta_ev = eps.stat.expected_value('theta')

# Get the mean value of positive theta
theta_evp = eps.stat.expected_value('theta',select={'theta':[0.0,None]})

# Get standard deviation of theta
theta_std = eps.stat.standard_deviation('theta')

# Get correlation coefficient between theta and ekin
theta_ekin_cc = eps.stat.correlation_coefficient('theta','ekin')

# Print informations
print('theta expected value : %.4E deg'%theta_ev)
print('positive theta expected value : %.4E deg'%theta_evp)
print('theta standard deviation : %.4E deg'%theta_std)
print('theta ekin correlation coefficient : %.4E'%theta_ekin_cc)
