# Particle Phase Space Analysis Toolkit

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

## Installation

Download the source code from github and add the following lines at the beginning of your script

```python
p2sat_path="/path/to/p2sat/"
import sys
if p2sat_path not in sys.path: sys.path.append(p2sat_path)

import p2sat
```

## Examples

Here is one quick example of p2sat usage, with my Geant4 app results (see ``examples/`` for more informations).

### Import results from a simulation file

```python
eps = p2sat.PhaseSpace(specie="electron")
eps.extract.txt("example.csv", sep=",")
```

### 1D histogram

Spectrum (Number/MeV) of all the electrons with time selection between 700 and 900 fs (bin width of 0.1 MeV)
```python
ekin,spectrum = eps.hist.h1('ekin',bwidth=0.1,select={'t':[700.0,900.0]})
```

### 1D histogram plot and fit

Spectrum of electrons, and exponential fit for energy > 0.511 MeV (bin width of 0.1 MeV, log scale)
```python
eps.plot.h1('ekin', log=True, bwidth=0.1)
eps.plot.f1('ekin', func_name="exp", log=True, bwidth=0.1, select={'ekin':[0.511,None]})
```

![](/mnt/local/esnault/Modules/p2sat/Figure_0.png)

### 2D histogram plot

Transverse particle dispersion (y, z) at x = 300 µm, for electrons with kinetic energy > 0.511 MeV (log color scale, bin width of 5 µm each, between -300 and 300 µm each)
```python
eps.plot.h2('y','z',log=True,
            bwidth1=5.0,bwidth2=5.0,
            brange1=[-300.,300.],brange2=[-300.,300.],
            select={'x':300,'ekin':[0.511,None]})
```

![](/mnt/local/esnault/Modules/p2sat/Figure_1.png)

