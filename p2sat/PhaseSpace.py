#coding:utf8
import numpy as np

from ._Data import _Data
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
  data : sub-object
    contains data data and methods to manipulate it, such as discretization or transformation.
  hist : sub-object
    make histograms from data data
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
      self.data     = _Data(self)
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

    w   = np.array([list(self.data.w) + list(other.data.w)])[0]
    x   = np.array([list(self.data.x) + list(other.data.x)])[0]
    y   = np.array([list(self.data.y) + list(other.data.y)])[0]
    z   = np.array([list(self.data.z) + list(other.data.z)])[0]
    px  = np.array([list(self.data.px) + list(other.data.px)])[0]
    py  = np.array([list(self.data.py) + list(other.data.py)])[0]
    pz  = np.array([list(self.data.pz) + list(other.data.pz)])[0]
    t   = np.array([list(self.data.t) + list(other.data.t)])[0]

    ps.data.update(w,x,y,z,px,py,pz,t)

    return ps

  def __len__(self):
    """
    Return the total number of configurations.
    """
    return len(self.data.w)


  def copy(self,verbose=False):
    """
    Return a copy of the current PhaseSpace object.

    Parameters
    ----------
    verbose : bool
      verbosity of the function. If True, a message is displayed when the attributes are loaded in memory
    """
    new = PhaseSpace()
    d=self.data
    new.data.update(d.w,d.x,d.y,d.z,d.px,d.py,d.pz,d.t,verbose=verbose)

    return new
