# Particle Phase Space Analysis Toolkit

p2sat is an open-source, object-oriented python package created to simplify particle phase-space analysis.

Basically, the particle phase-space informations (statistical weight, position, momentum ; one line for one particle) are stored in numpy arrays, and some methods are available to make 1d/2d/3d histograms and plots in a simple way. Other particles informations (kinetic energy, divergence angles, ...) are automatically calculated from particle phase-space and any of them can be used to filter an histogram or a plot, such as filtering a divergence plot to represent only particles that are in a defined position, or having a defined energy range.

See objects documentation for more informations.

**Note :**
This package was made for my personal use and then contains only few methods to import data from code results, but thank to object inheritance, it can be adapted very easily to other codes (see `PhaseSpaceGeneric` object in *p2sat/PhaseSpace.py*).

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

```python
import matplotlib.pyplot as plt

es = p2sat.PhaseSpaceGeant4()
es.extract("../Geant4/ps_nt_electron_t*.csv",nthreads=10)

es.plot.h1('ekin',bwidth=0.1,select={'x':[50.0,100.0]})
plt.title('Electron spectrum for $x \in [50,100] \mu m$')
plt.xlabel('ekin (MeV)')
plt.ylabel('Number per MeV')

es.plot.h2('y','z',bwidth1=10.0,bwidth2=10.0,brange1=[-500.,500.],brange2=[-500.,500.],select={'x':50,'ekin':[0.511,None]})
plt.title('Transverse particle dispersion at $x=50$, for $E_{kin} > 0.511 MeV$')
plt.xlabel('y (um)')
plt.ylabel('z (um)')
```
