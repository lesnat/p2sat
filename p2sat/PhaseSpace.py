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
    contains raw data and methods to manipulate it, such as discretization or transformation.
  hist : sub-object
    make histograms from data
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
      if specie in ("gamma","g"):
          self.specie = "gamma"
          self.mass = 0
      elif specie in ("positron","e+"):
          self.specie = "e+"
          self.mass = 0.511
      elif specie in ("electron","e-"):
          self.specie = "e-"
          self.mass = 0.511
      elif specie in ("muon+","mu+"):
          self.specie = "mu+"
          self.mass = 105.6
      elif specie in ("muon-","mu-"):
          self.specie = "mu-"
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

  def __str__(self):
    txt  = "\n"
    txt += "p2sat PhaseSpace instance located at %s\n\n"%hex(id(self))
    txt += "Specie                      : %s\n"%self.specie
    txt += "Number of configurations    : %i\n"%len(self)
    txt += "Total number of particles   : %.4E\n\n"%sum(self.data.w)
    txt += "Statistics  : ( min      ,  max      ,  mean     ,  std      ) unit\n"
    txt += "    w       : (% .3E, % .3E, % .3E, % .3E) \n"%(min(self.data.w),max(self.data.w),self.data.w.mean(),self.data.w.std())
    txt += "    t       : (% .3E, % .3E, % .3E, % .3E) fs\n\n"%(min(self.data.t),max(self.data.t),self.data.t.mean(),self.data.t.std())

    txt += "    x       : (% .3E, % .3E, % .3E, % .3E) um\n"%(min(self.data.x),max(self.data.x),self.data.x.mean(),self.data.x.std())
    txt += "    y       : (% .3E, % .3E, % .3E, % .3E) um\n"%(min(self.data.y),max(self.data.y),self.data.y.mean(),self.data.y.std())
    txt += "    z       : (% .3E, % .3E, % .3E, % .3E) um\n"%(min(self.data.z),max(self.data.z),self.data.z.mean(),self.data.z.std())
    txt += "    r       : (% .3E, % .3E, % .3E, % .3E) um\n\n"%(min(self.data.r),max(self.data.r),self.data.r.mean(),self.data.r.std())

    txt += "    px      : (% .3E, % .3E, % .3E, % .3E) MeV/c\n"%(min(self.data.px),max(self.data.px),self.data.px.mean(),self.data.px.std())
    txt += "    py      : (% .3E, % .3E, % .3E, % .3E) MeV/c\n"%(min(self.data.py),max(self.data.py),self.data.py.mean(),self.data.py.std())
    txt += "    pz      : (% .3E, % .3E, % .3E, % .3E) MeV/c\n"%(min(self.data.pz),max(self.data.pz),self.data.pz.mean(),self.data.pz.std())
    txt += "    p       : (% .3E, % .3E, % .3E, % .3E) MeV/c\n\n"%(min(self.data.p),max(self.data.p),self.data.p.mean(),self.data.p.std())

    txt += "    ekin    : (% .3E, % .3E, % .3E, % .3E) MeV\n"%(min(self.data.ekin),max(self.data.ekin),self.data.ekin.mean(),self.data.ekin.std())
    txt += "    gamma   : (% .3E, % .3E, % .3E, % .3E) \n"%(min(self.data.gamma),max(self.data.gamma),self.data.gamma.mean(),self.data.gamma.std())
    txt += "    beta    : (% .3E, % .3E, % .3E, % .3E) \n"%(min(self.data.beta),max(self.data.beta),self.data.beta.mean(),self.data.beta.std())
    txt += "    v       : (% .3E, % .3E, % .3E, % .3E) um/fs\n\n"%(min(self.data.v),max(self.data.v),self.data.v.mean(),self.data.v.std())

    txt += "    theta   : (% .3E, % .3E, % .3E, % .3E) deg\n"%(min(self.data.theta),max(self.data.theta),self.data.theta.mean(),self.data.theta.std())
    txt += "    phi     : (% .3E, % .3E, % .3E, % .3E) deg\n"%(min(self.data.phi),max(self.data.phi),self.data.phi.mean(),self.data.phi.std())
    txt += ""

    return txt

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
