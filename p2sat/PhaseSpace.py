#coding:utf8
import numpy as np

from ._Raw import _Raw
from ._Hist import _Hist
from ._Plot import _Plot
from ._Export import _Export
from ._Extract import _Extract
from ._Stat import _Stat


class PhaseSpace(object):
  """
  Base class for particle phase-space analysis.

  Parameters
  ----------
  specie : str
    Name of the particle specie. Availables are gamma,e-,e+,mu-,mu+.

  Attributes
  ----------
  raw : sub-object
    contains raw data and methods to manipulate it, such as discretization or transformation.
  hist : sub-object
    make histograms from raw data
  plot : sub-object
    plot histograms
  extract : sub-object
    load phase space from a file
  export : sub-object
    export phase space into a file
  stat : sub-object
    make statistics on particle phase space

  Notes
  -----
  See sub-objects documentation for more informations
  """
  def __init__(self,specie):
      self.specie = specie
      if specie in ("gamma","g"):
          self.mass = 0
      elif specie in ("positron","electron","e-","e+","e"):
          self.mass = 0.511
      elif specie in ("muon","muon+","muon-","mu+","mu-","mu"):
          self.mass = 105.6
      else:
          raise NameError("Unknown particle specie.")
      self.raw      = _Raw(self)
      self.hist     = _Hist(self)
      self.plot     = _Plot(self)
      self.extract  = _Extract(self)
      self.export   = _Export(self)
      self.stat     = _Stat(self)

  def __add__(self,other):
    """
    Return a new PhaseSpace object, combination of the 2 previous.
    """
    ps = PhaseSpace()

    w   = np.array([list(self.raw.w) + list(other.raw.w)])[0]
    x   = np.array([list(self.raw.x) + list(other.raw.x)])[0]
    y   = np.array([list(self.raw.y) + list(other.raw.y)])[0]
    z   = np.array([list(self.raw.z) + list(other.raw.z)])[0]
    px  = np.array([list(self.raw.px) + list(other.raw.px)])[0]
    py  = np.array([list(self.raw.py) + list(other.raw.py)])[0]
    pz  = np.array([list(self.raw.pz) + list(other.raw.pz)])[0]
    t   = np.array([list(self.raw.t) + list(other.raw.t)])[0]

    ps.raw.update(w,x,y,z,px,py,pz,t)

    return ps

  def __len__(self):
    """
    Return the total number of configurations.
    """
    return len(self.raw.w)


  def copy(self,verbose=False):
    """
    Return a copy of the current PhaseSpace object.

    Parameters
    ----------
    verbose : bool
      verbosity of the function. If True, a message is displayed when the attributes are loaded in memory
    """
    new = PhaseSpace()
    r=self.raw
    new.raw.update(r.w,r.x,r.y,r.z,r.px,r.py,r.pz,r.t,verbose=verbose)

    return new
