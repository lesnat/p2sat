============
Installation
============

The most simple way to install p2sat is to use pip (https://pypi.org/project/p2sat/)
::
  pip install p2sat

Otherwise, you can also download the source code from github (https://github.com/lesnat/p2sat), extract it and and type the following commands
::
  cd p2sat
  python setup.py install

If it is not working, you can add the following lines at the beginning of your script
::
  p2sat_path="/path/to/p2sat/"
  import sys
  if p2sat_path not in sys.path: sys.path.append(p2sat_path)

  import p2sat


p2sat is written for python 2.7 but might be compatible with 3.

Its only dependancies are python packages numpy and matplotlib.
