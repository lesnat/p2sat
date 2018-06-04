#coding:utf8
import numpy as np

class _Raw(object):
  """
  Class containing raw data and methods to manipulate it.

  Attributes
  ----------
  w : numpy.ndarray
    particle statistical weight
  x,y,z : numpy.ndarray
    particle x,y,z position in um
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
  - ekin is defined as :math:`(\sqrt{(p/m_e c)^2+1}-1) \\times m_e c^2` *
  - theta is defined as :math:`\\arctan{p_y/p_x}`
  - phi is defined (yet) as :math:`\\arctan{p_z/p_x}`


  * detail of the calculus can be found at ... TODO
  """
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
    self.ekin   = (np.sqrt((self.p/0.511)**2 + 1) - 1) * 0.511
    self.gamma  = self.ekin/0.511 + 1.

    """
    self.w.flags.writeable  = False

    self.x.flags.writeable  = False
    self.y.flags.writeable  = False
    self.z.flags.writeable  = False

    self.px.flags.writeable = False
    self.py.flags.writeable = False
    self.pz.flags.writeable = False

    self.r.flags.writeable  = False
    self.p.flags.writeable  = False
    self.theta.flags.writeable = False
    self.phi.flags.writeable = False
    self.ekin.flags.writeable = False
    """
    if verbose: print("Done !")
