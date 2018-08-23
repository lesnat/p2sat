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
es = p2sat.PhaseSpace(specie="electron")
es.extract.Geant4_csv("Al_target",nthreads=10)
```

### 1D histogram

Spectrum (Number/MeV) of all the electrons with time selection between 0 and 100 fs (bin width of 0.1 MeV)
```python
ekin,spectrum = es.hist.h1('ekin',bwidth=0.1,select={'t':[0.0,100.0]})
```

### 1D histogram plot and fit

Spectrum of electrons with x position between 50 and 100 µm (bin width of 0.1 MeV, exponential fit, log scale)
```python
es.plot.h1('ekin',log=True,bwidth=0.1,select={'x':[50.0,100.0]})
es.plot.f1('ekin',func_name="exp",log=True,bwidth=0.1,select={'x':[50.0,100.0]})
```

### 2D histogram plot

Transverse particle dispersion at x = 50 µm, for electrons with kinetic energy > 0.511 MeV (bin width of 10 µm between -500 and 500 µm)
```python
es.plot.h2('y','z',bwidth1=10.0,bwidth2=10.0,brange1=[-500.,500.],brange2=[-500.,500.],select={'x':50.0,'ekin':[0.511,None]})
```
