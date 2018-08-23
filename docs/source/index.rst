.. p2sat documentation master file, created by
   sphinx-quickstart on Thu Aug 23 14:20:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


*****
p2sat
*****

p2sat is an open-source, object-oriented python package created to simplify particle phase-space analysis.

Core features of the package are :

- Automatic calculation of kinetic energy, divergence angle and gamma factor of the particles from phase space informations
- Histogram making (1D, 2D, 3D) and data fits (1D)
- Plotting (1D to 3D histograms, scatter and contour plots) with automatic normalizations and legend
- Particle filtering with a given property (for example select all the particles at a given position)
- Statistical tools (standard deviation, covariance, ...)
- Import data from simulation files (Smilei, Geant4, text files, ...)
- Low memory load

This allows to process complex operations in a very concise and clear way, as shown in the examples.

See objects documentation for more informations.

**Notes :**

- This package was made for my personal use and then contains only few methods to import data from code results, but you can easily add your own (please, share !) and use it to perform your data analysis. See sub-object ``_Extract`` for more informations.
- This tool can be usefull to physicists working with Particle-In-Cell or Monte Carlo codes

Contents
--------

.. toctree::
  :maxdepth: 3

  intro
  install
  use
  examples
