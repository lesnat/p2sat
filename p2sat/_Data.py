#coding:utf8
import numpy as np

class _Data(object):
  """
  Class containing raw data and methods to manipulate it.

  Attributes
  ----------
  w : numpy.ndarray
    particle statistical weight
  x,y,z : numpy.ndarray
    particle x,y,z position in um
  t : numpy.ndarray
    particle time in fs
  r : numpy.ndarray
    absolute distance to the x axis in um
  px,py,pz : numpy.ndarray
    particle momentum in x,y,z direction in MeV/c
  p : numpy.ndarray
    absolute particle momentum in MeV/c
  ekin : numpy.ndarray
    particle energy in MeV
  theta : numpy.ndarray
    angle between px and py in degree
  phi : numpy.ndarray
    angle between ??? in degree

  Notes
  -----
  As all the calculations are done with the previously defined units,
  the input data might be firstly converted to those units.

  Calculations :

  - r is defined as :math:`\sqrt{y^2+z^2}`
  - p is defined as :math:`\sqrt{p_x^2+p_y^2+p_z^2}`
  - theta is defined as :math:`\\arctan{p_y/p_x}`
  - phi is defined (yet) as :math:`\\arctan{p_z/p_x}`
  - ekin is defined as

    - :math:`(\sqrt{(p/m_e c)^2+1}-1) \\times m_e c^2` for massive species
    - :math:`p` otherwise (here ekin is the total particle energy)
  - gamma is defined as

    - :math:`E_{kin}/m_e c^2 + 1` for massive species
    - :math:`...` otherwise

  Details of the calculations can be found at ... TODO
  """
  def __init__(self,PhaseSpace):
    self._ps= PhaseSpace
    self.w  = None
    self.x  = None
    self.y  = None
    self.z  = None
    self.px = None
    self.py = None
    self.pz = None
    self.t  = None

  def update(self,w,x,y,z,px,py,pz,t,verbose=True):
    """
    Update class attributes with new values.

    Parameters
    ----------
    w,x,y,z,px,py,pz : list or numpy.ndarray
      particle phase space. More information can be found in raw object documentation
    verbose : bool
      verbosity of the function. If True, a message is displayed when the attributes are loaded in memory

    TODO: get np array to be immutable with x.writeable=False ?
    """
    if verbose: print("Updating raw values ...")
    # Save values into np.array objects
    self.w  = np.array(w)
    self.x  = np.array(x)
    self.y  = np.array(y)
    self.z  = np.array(z)
    self.px = np.array(px)
    self.py = np.array(py)
    self.pz = np.array(pz)
    self.t  = np.array(t)

    # Calculate other parameters from it
    self.r      = np.sqrt(self.y**2+self.z**2)
    self.p      = np.sqrt(self.px**2+self.py**2+self.pz**2)
    self.theta  = np.degrees(np.arctan(self.py/self.px))
    self.phi    = np.degrees(np.arctan(self.pz/self.px)) # Geometrical effect ? change -> 0 pi
    mass = self._ps.mass
    if mass == 0:
        self.ekin   = self.p
        self.gamma  = self.ekin # FIXME : vérifier
    else:
        self.ekin   = (np.sqrt((self.p/mass)**2 + 1) - 1) * mass
        self.gamma  = self.ekin/mass + 1.
    if verbose: print("Done !")

  def generate(self,**kargs):
      """
      Generate a particle phase space from given laws

      TODO
      """
      pass

  def select(self,axis,faxis,frange,fpp=1e-7):
    """
    Filter an axis with a value/range on another axis.

    Parameters
    ----------
    axis : str or numpy.ndarray
      axis to filter
    faxis : list of str or list of numpy.ndarray
      filtering axis
    frange : list of int, float, list/tuple of 2 float
      filtering value/range (value if int, range if float or list/tuple). If a frange element is None, the minimum/maximum value is taken
    fpp : float, optional
      relative floating point precision. Default is 1e-7

    Returns
    -------
    axis : numpy.ndarray
      filtered axis

    Examples
    --------

    It is possible to filter by an int value

    >>> w = np.random.uniform(low=0.,high=10.,size=10)
    >>> x = np.array([1,3,3,3,7,9,5,3,7,3])
    >>> w = select(w,x,3) # Select all the w satisfying x==3

    or filter by a range

    >>> w = np.random.uniform(low=0.,high=10.,size=1000)
    >>> ekin = np.random.exponential(scale=3.0,size=1000)
    >>> w = select(w,ekin,[0.511,None]) # Select all the w with :math:`ekin \in [0.511,+\infty] MeV`

    If frange is a list/tuple or a float, the filtering is done with a fpp precision
    """
    # if type(faxis) is not (list or tuple): faxis=[faxis]

    # Get a copy of axis and faxis (from a str or not)
    if type(axis) is str: axis=eval("self.%s"%axis)
    axis=np.array(axis)
    for i,fax in enumerate(faxis):
      if type(fax) is str: faxis[i]=eval("self.%s"%faxis[i])
      faxis[i]=np.array(faxis[i])

    # Filtering ...
    for i,_ in enumerate(faxis):
      if type(frange[i]) is list or type(frange[i]) is tuple:
        if frange[i][0] is None: frange[i][0]=min(axis)
        if frange[i][1] is None: frange[i][1]=max(axis)
        filtr=np.array([x>frange[i][0]*(1-fpp) and x<frange[i][1]*(1+fpp) for x in faxis[i]])
        axis=axis[filtr]
      elif type(frange[i]) is int:
        axis=axis[faxis[i]==frange[i]]
      elif type(frange[i]) is float:
        axis=self.select(axis,faxis=[faxis[i]],frange=[[frange[i],frange[i]]])
      else:
        raise TypeError('frange type must be int/float or list/tuple of 2 float.')

      # filter next faxis with current faxis
      if len(faxis)>i+1:faxis[i+1]=self.select(faxis[i+1],[faxis[i]],[frange[i]])

    return axis

  def transformate(self,translation=(0.0,0.0,0.0),rotation=(0.0,0.0)):
    """
    Transformate the particle phase space with given translation and rotation.

    TODO
    """
    pass

  def propagate(self):
    """
    Propagate the phase space to a given position or time.

    TODO
    """
    pass

  def discretize(self,with_time=True,verbose=True,**kargs):
    """
    Discretize the particles phase space in a 6 or 7 D histogram.

    Parameters
    ----------
    with_time : bool, optional
      discretize with time (7D). Default is True
    verbose : bool, optional
      verbosity. Default is True
    kargs
      optional keyword arguments to pass to the hist.hn function

    Notes
    -----
    This method can be used to significantly reduce disk space usage
    when saving data into output file.

    See Also
    --------
    hn
    """
    hn=self._ps.hist.hn

    if verbose : print('Data discretization ...')
    #bi,hi=hn(['x','y','z','px','py','pz'],bwidth=bwidth,brange=brange,wnorm=[1.0]*6,select=select)
    if with_time:
      bi,hi=hn([self.x,self.y,self.z,self.px,self.py,self.pz,self.t],wnorm=[1.0]*7,**kargs)
      bx,by,bz,bpx,bpy,bpz,bt=bi
    else:
      bi,hi=hn([self.x,self.y,self.z,self.px,self.py,self.pz],wnorm=[1.0]*6,**kargs)
      bx,by,bz,bpx,bpy,bpz=bi
    if verbose : print('Done !')
    w       = []
    x,y,z   = [],[],[]
    px,py,pz= [],[],[]
    t       = []

    if verbose : print('Getting new configurations ...')
    hb=hi.nonzero()

    print(bx,by,bz,bpx,bpy,bpz)

    w   = hi[hb]
    x   = bx[hb[0]]
    y   = by[hb[1]]
    z   = bz[hb[2]]
    px  = bpx[hb[3]]
    py  = bpy[hb[4]]
    pz  = bpz[hb[5]]
    if with_time:
      t   = bt[hb[6]]
    else:
      t   = [0.0]*len(w)

    if verbose : print('Done !')

    self.update(w,x,y,z,px,py,pz,t,verbose=verbose)
