============
Installation
============

p2sat is written for python 2.7 but might be compatible with 3.

Its only dependancies are python packages `numpy` and `matplotlib`.

Download and extract the source code, and use the following lines at the beginning of your script

```python
import sys
p2sat_path='/path/to/p2sat/'
if p2sat_path not in sys.path:sys.path.append(p2sat_path)
import p2sat
```
