#coding:utf8
import numpy as np

from ._Raw import _Raw
from ._Hist import _Hist
from ._Plot import _Plot
from ._Tools import _Tools
from ._Extract import _Extract


class PhaseSpace(object):
  """
  Base class for particle phase-space analysis.
  
  Attributes
  ----------
  raw : sub-object
    contains raw data and methods to manipulate it, such as import/export into a file or value filtering method. Appropriate place to IO ?
  hist : sub-object
    contains methods to make histograms from raw data
  plot : sub-object
    contains methods to plot histos
  extract : sub-object
    load phase space from a simulation file
  
  Notes
  -----
  See sub-objects documentation for more informations
  """
  def __init__(self):
      self.raw  = _Raw()
      self.hist = _Hist(self)
      self.plot = _Plot(self)
      self.tools= _Tools(self)
      self.extract=_Extract(self)    





######### DEPRECATED ############

class PhaseSpaceGeneric(object):
  """
  Base class for particle phase-space analysis.
  
  Attributes
  ----------
  extract : "polymorphic" method
    load phase space from a simulation file
  raw : sub-object
    contains raw data and methods to manipulate it, such as import/export into a file or value filtering method. Appropriate place to IO ?
  hist : sub-object
    contains methods to make histograms from raw data
  plot : sub-object
    contains methods to plot histos
  
  Notes
  -----
  See sub-objects documentation for more informations
  """
  def __init__(self):
      self.raw  = _Raw()
      self.hist = _Hist(self)
      self.plot = _Plot(self)
      self.tools= _Tools(self)
  
  def extract(self):
    """
    Extract raw data from a simulation file.
    
    Notes
    -----
    This abstract method must be overwritten by the PhaseSpaceGeneric child classes
    """
    raise NotImplementedError
      

class PhaseSpaceSmilei(PhaseSpaceGeneric):
  def __init__(self):
      self.raw  = _Raw()
      self.hist = _Hist(self)
      self.plot = _Plot(self)
      self.tools= _Tools(self)
      self.__doc__=PhaseSpaceGeneric.__doc__
  
  def extract(self,Screen,timestep,xnorm=1.0,wnorm=1.0,X=0.0,verbose=True):
    """
    Extract simulation results from a Smilei Screen diagnostic
    
    Parameters
    ----------
    Screen : Smilei object
      simulation diagnostic to open
    timestep : int
      timestep at which open the diag
    xnorm : float
      normalization factor of x axis (laser propagation axis), into um
    wnorm : float
      normalization factor of weight, into physical number of particles
    X : float
      Screen diagnostic position from the beginning of target solid density, in um
    verbose : bool, optional
      verbosity
    
    Notes
    -----
    The Smilei Screen diagnostic should be declared in the namelist at the following format
      1D:
      DiagScreen( ...,
            axis = [
                  ['py'   , -pmax/10, pmax/10 , 301],
                  ['px'   , pmin    , pmax    , 301]
                  ]
      )
      
      2D:
      DiagScreen( ...,
            axis = [
                  ['y'    , -rmax   , rmax    , 301],
                  ['pz'   , -pmax/10, pmax/10 , 301],
                  ['py'   , -pmax/10, pmax/10 , 301],
                  ['px'   , pmin    , pmax    , 301]
                  ]
      )
      
      3D :
      DiagScreen( ...,
            axis = [
                  ['z'    , -rmax   , rmax    , 301],
                  ['y'    , -rmax   , rmax    , 301],
                  ['pz'   , -pmax/10, pmax/10 , 301],
                  ['py'   , -pmax/10, pmax/10 , 301],
                  ['px'   , pmin    , pmax    , 301]
                  ]
      )
      
    The order is very important. It is also recommended to have an odd number of bins.
    
    The weights normalization only works for linearly spaced bins.
    
    Examples
    --------
    >>> import happi
    >>> S = happi.Open('../Smilei/')
    >>> nl = S.namelist
    >>> es = p2sat.PhaseSpaceSmilei()
    >>> es.extract(S.Screen.Screen2,nl.Ndt,xnorm=nl.um,wnorm=nc*vol,X=1.0) # With ``nc`` the critical density & ``vol`` the characteristic interaction volume
    
    TODO
    ----
    - Implement
    - Check if the format is OK
    - Check if giving X is mandatory (another way to do it with nl.Screen ?)
    """
    if verbose:print("Extracting screen data ...")
    
    w         = []
    x,y,z     = [],[],[]
    px,py,pz  = [],[],[]
    wNorm     = float(wnorm) # Get a copy of wnorm
    
    dpx,dpy,dpz=[],[],[]
    
    data  = Screen(timesteps=timestep).getData()[0] # TODO: check this
    
    try:
      Z     = Screen().get()['z']
      wNorm*= (max(Z)-min(Z))/len(Z)
      Z    *= xnorm
    except KeyError:
      Z     = [0.0]
      data  = [data]
    
    try:
      Y     = Screen().get()['y']
      wNorm*= (max(Y)-min(Y))/len(Y)
      Y    *= xnorm
    except KeyError:
      Y     = [0.0]
      data  = [data]
    
    try:
      Pz    = Screen().get()['pz']
      wNorm*= (max(Pz)-min(Pz))/len(Pz)
      Pz   *= 0.511
    except KeyError:
      Pz    = [0.0]
      data  = [data]
    
    Py      = Screen().get()['py']
    Px      = Screen().get()['px']
    wNorm  *= (max(Py)-min(Py))/len(Py)
    wNorm  *= (max(Px)-min(Px))/len(Px)
    Py     *= 0.511
    Px     *= 0.511
    
    for iz,ez in enumerate(data):
      for iy,ey in enumerate(ez):
        for ipz,epz in enumerate(ey):
          for ipy,epy in enumerate(epz):
            for ipx,epx in enumerate(epy):
              if epx!=0:
                x.append(X)
                y.append(Y[iy])
                z.append(Z[iz])
                px.append(Px[ipx])
                py.append(Py[ipy])
                pz.append(Pz[ipz])
                w.append(epx*wNorm)
                dpx.append(Px[ipx+1] - Px[ipx])
                dpy.append(Py[ipy+1] - Py[ipy])
                if len(Pz)==1:
                  dpz.append(0.0)
                else:
                  dpz.append(Pz[ipz+1] - Pz[ipz])
              
    if verbose:print("Done !")

    self.raw.update(w,x,y,z,px,py,pz,verbose)
    
    self.raw.dpx=np.array(dpx)
    self.raw.dpy=np.array(dpy)
    self.raw.dpz=np.array(dpz)
    
    #self.raw.dp = np.sqrt(self.raw.dpx**2 + self.raw.dpy**2 + self.raw.dpz**2)
    #self.raw.dekin = (np.sqrt((self.raw.dp/0.511)**2+1.0)-1.0)*0.511
    self.raw.dp     = self.raw.px/self.raw.p * self.raw.dpx + self.raw.py/self.raw.p * self.raw.dpy + self.raw.pz/self.raw.p * self.raw.dpz
    self.raw.dekin  = self.raw.p/np.sqrt(self.raw.p**2 + 1.0) * self.raw.dp
  
  def extract_1d(self,Screen,timestep,xnorm,wnorm=1.0,X=0):
    w         = []
    x,y,z     = [],[],[]
    px,py,pz  = [],[],[]
    
    data= Screen(timesteps=timestep).getData()[0]
    Px  = Screen().get()['px'] * 0.511
    Py  = Screen().get()['py'] * 0.511
    
    try:
      Y = Screen().get()['y'] * xnorm
    except:
      Y = [0.0]
      data = [data]
    
    wNorm = wnorm/(0.511**2)
    wNorm *= (max(Px)-min(Px))/len(Px)
    wNorm *= (max(Py)-min(Py))/len(Py)
    
    # redéfinir l'ordre dans namelist ?
    print("Extracting screen data ...")
    for iy,ey in enumerate(data):
      for ipx,epx in enumerate(ey):
        for ipy,epy in enumerate(epx):
          if epy!=0:
            y.append(Y[iy])
            px.append(Px[ipx])
            py.append(Py[ipy])
            w.append(epy*wNorm)

    pz = [0.0] * len(w)
    x = [X] * len(w)
    z = [0.0] * len(w)
    print("Data succesfully imported")
    self.raw.update(w,x,y,z,px,py,pz)
    
  def extract_2d(self,Screen,timestep,xnorm,wnorm=1.0,X=0):
    w         = []
    x,y,z     = [],[],[]
    px,py,pz  = [],[],[]
    
    data= Screen(timesteps=timestep).getData()[0]
    Px  = Screen().get()['px'] * 0.511
    Py  = Screen().get()['py'] * 0.511
    Pz  = Screen().get()['pz'] * 0.511


    
    try:
      Y = Screen().get()['y'] * xnorm
    except:
      Y = [0.0]
      data = [data]
    
    wNorm  = wnorm/(0.511**3)
    wNorm *= (max(Px)-min(Px))/len(Px)
    wNorm *= (max(Py)-min(Py))/len(Py)
    wNorm *= (max(Pz)-min(Pz))/len(Pz)
    # redéfinir l'ordre dans namelist ?
    print("Extracting screen data ...")
    for iy,ey in enumerate(data):
      for ipx,epx in enumerate(ey):
        for ipy,epy in enumerate(epx):
          for ipz,epz in enumerate(epy):
            if epz!=0:
              y.append(Y[iy])
              px.append(Px[ipx])
              py.append(Py[ipy])
              pz.append(Pz[ipz])
              w.append(epz*wNorm)
              
    x = [X] * len(w)
    z = [0.0] * len(w)
    print("Data succesfully imported")
    self.raw.update(w,x,y,z,px,py,pz)


class PhaseSpaceGeant4(PhaseSpaceGeneric):
  def __init__(self):
    self.raw  = _Raw()
    self.hist = _Hist(self)
    self.plot = _Plot(self)
    self.tools= _Tools(self)
    self.__doc__=PhaseSpaceGeneric.__doc__
    
  def extract(self,file_name,nthreads=1,verbose=True):
    """
    Extract simulation results from a Geant4 NTuple output file
    
    Parameters
    ----------
    file_name : str
      name of the output file. If it ends with '*_t0.*', the number '0' will be replaced by the number of the current thread
    nthreads : int
      total number of threads to consider
    verbose : bool, optional
      verbosity
    
    Notes
    -----
    The Geant4 NTuple format should be
      w,x,y,z,px,py,pz
      . . . . .  .  .
      . . . . .  .  .
      . . . . .  .  .
      
    Yet only xml and csv file format are accepted
    
    Examples
    --------
    >>> eg = p2sat.PhaseSpaceGeant4()
    >>> eg.extract("../Geant4/testem_nt_electron_t*.csv",nthreads=10)
    
    TODO
    ----
    - while True + try/except to loop over nthreads ?
    """
    data = []
    fext = file_name.split('.')[-1]   # File extension
    fbase= file_name[:-(len(fext)+1)] # File base name
    
    for thread in range(0,nthreads):
      fname=fbase[:-1]+str(thread)+"."+fext
      if verbose:print("Extracting %s ..."%fname)

      if fext=="csv":
        with open(fname,'r') as f:
          for line in f.readlines():
            if line[0]!='#':
              for e in line.split(','):
                data.append(float(e))
      elif fext=="xml":
        from lxml import etree      
        with etree.parse(fname) as tree:
          for entry in tree.xpath('/aida/tuple/rows/row/entry'):
            data.append(float(entry.get('value')))
      else:
        raise NameError("Unknown file extension : %s"%fext)
    
    w   = data[0::7]
    x   = data[1::7]
    y   = data[2::7]
    z   = data[3::7]
    px  = data[4::7]
    py  = data[5::7]
    pz  = data[6::7]
    if verbose:print("Done !")
    
    self.raw.update(w,x,y,z,px,py,pz,verbose)


class PhaseSpaceTrILEns(PhaseSpaceGeneric):
  def __init__(self):
    self.raw  = _Raw()
    self.hist = _Hist(self)
    self.plot = _Plot(self)
    self.tools= _Tools(self)
    self.__doc__=PhaseSpaceGeneric.__doc__
    
  def extract(self,path,specie,verbose=True):
    """
    Extract simulation results from a TrILEns output.txt file
    
    Parameters
    ----------
    path : str
      simulation path
    specie : str
      specie to find in the output. The specie name must be in plural form (i.e 'electrons' or 'positrons')
    verbose : bool, optional
      verbosity
    
    Notes
    -----
    ...
    
    Examples
    --------
    >>> et = p2sat.PhaseSpaceTrILEns()
    >>> et.extract("../TrILEns/",specie="positrons")
    
    TODO
    ----
    - Check Interface.f90 at line ~ 120 -> condition for e-/e+ & not working with photons
    """
    if verbose:print("Extracting {} data from {}output.txt ...".format(specie,path))
    
    w         = []
    x,y,z     = [],[],[]
    px,py,pz  = [],[],[]
    
    is_correct_specie=False
    
    with open(path+'output.txt','r') as f:
      for _ in range(34):
        f.readline()
      
      for line in f.readlines():
        try:
          W,X,Y,Z,Px,Py,Pz,Gamma,Chi=line.split()
          if is_correct_specie:
            w.append(float(W))
            x.append(float(X))     ; y.append(float(Y))   ; z.append(float(Z))
            px.append(float(Px))   ; py.append(float(Py)) ; pz.append(float(Pz))
        except ValueError:
          if specie in line.split():
            is_correct_specie = True
          else:
            is_correct_specie = False
    
    if verbose:print("Data succesfully imported")
    
    self.raw.update(w,x,y,z,px,py,pz,verbose)
