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

This allows to plot complicated graphs in a very concise and clear way, as shown in the following examples.

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

## Example

Here is one quick example of p2sat plot usage, with my Geant4 app results.

### Import results from a simulation file
```python
# Instanciate a PhaseSpace object for electron
es = p2sat.PhaseSpace(specie="electron")
# Import results from a simulation file
es.extract.Geant4_csv("Al_target",nthreads=10)
```

### Plot and fit electron spectrum for x position between 50 and 100 µm
```python
es.plot.h1('ekin',bwidth=0.1,select={'x':[50.0,100.0]})
es.plot.f1('ekin',bwidth=0.1,select={'x':[50.0,100.0]})
```

### Plot transverse particle dispersion at x = 50µm, for electrons with energy > 0.511 MeV
```python
es.plot.h2('y','z',bwidth1=10.0,bwidth2=10.0,brange1=[-500.,500.],brange2=[-500.,500.],select={'x':50,'ekin':[0.511,None]})
```
