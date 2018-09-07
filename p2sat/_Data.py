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
    - :math:`+\infty` otherwise

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
    self.labels = {
              'x'     : '$x$',
              'y'     : '$y$',
              'z'     : '$z$',
              'px'    : '$p_x$',
              'py'    : '$p_y$',
              'pz'    : '$p_z$',
              't'     : '$t$',

              'r'     : '$r$',
              'p'     : '$p$',
              'ekin'  : '$E_{kin}$',
              'gamma' : '$\gamma$',

              'theta' : '$\\theta$',
              'phi'   : '$\\phi$'
             }

    self.units = {
              'x'     : 'um',
              'y'     : 'um',
              'z'     : 'um',
              'px'    : 'MeV/c',
              'py'    : 'MeV/c',
              'pz'    : 'MeV/c',
              't'     : 'fs',

              'r'     : 'um',
              'p'     : 'MeV/c',
              'ekin'  : 'MeV',
              'gamma' : None,

              'theta' : 'deg',
              'phi'   : 'deg'
              }

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
        self.gamma  = np.inf
    else:
        self.ekin   = (np.sqrt((self.p/mass)**2 + 1) - 1) * mass
        self.gamma  = self.ekin/mass + 1.
    if verbose: print("Done !")

  def generate(self,Nconf,Npart,ekin,theta,phi,pos=None,time=None,verbose=True):
      """
      Generate a particle phase space from given laws.

      Parameters
      ----------
      Nconf : int
        total number of configurations
      Npart : float
        total number of particles
      ekin : dict
        parameters to generate kinetic energy
      theta : dict
        parameters to generate theta angle distribution
      phi : dict
        parameters to generate phi angle distribution
      pos : dict, optional
        parameters to generate position distribution. Default is 0 for x,y,z
      time : dict
        parameters to generate time distribution. Default is 0

      Notes
      -----
      The dictionnaries must each time at least contain the key 'law' with a value
      depending on which law are available for each physical quantity

      For dict `ekin`, available laws are :
      - 'mono', for a mono-energetic source. Energy must be given as a value of keyword 'energy'
      - 'exp', for exponential energy. Temperature must be given as a value of keyword 'T'

      For dict `theta` and `phi`, available laws are :
      - 'mono', for a directional source. Angle must be given as a value of keyword 'angle'
      - 'iso', for an isotropic source. An optional keyword 'max' can be given to specify a maximum angle
      - 'gauss', for a gaussian spreading. Center of the distribution must be given with keyword 'mu', and standard deviantion with keyword 'sigma'

      Examples
      --------
      Assuming a `PhaseSpace` object is instanciated as `eps`, you can generate
      a mono-energetic source in isotropic direction for 1e12 particles represented
      by 1e6 configurations as follows

      >>> eps.data.generate(Nconf=1e6,Npart=1e12,ekin={"law":"mono","E":20.0},theta={"law":"iso"},phi={"law":"iso"})
      """
      # Print a starting message
      if verbose:
        print("Generate particle phase-space for \"%s\" ekin law, \"%s\" theta law, \"%s\" phi law ..."
        %(ekin["law"],theta["law"],phi["law"]))

      # Ensure that Nconf is of type int (for values such as 1e6)
      Nconf = int(Nconf)
      # Generate weights
      wnorm = float(Npart)/Nconf
      w = np.array([wnorm] * Nconf)

      # Generate theta angle
      if theta["law"]=="mono":
        theta0 = np.array([theta["angle"]] * Nconf)
      elif theta["law"]=="iso":
        try:
          mangle = theta["max"]*np.pi/180.
        except KeyError:
          mangle = np.pi
        theta0 = np.random.uniform(0.,mangle,Nconf)
      elif theta["law"]=="gauss":
        theta0 = np.random.normal(theta["mu"],theta["sigma"],Nconf)
      # Generate phi angle
      if phi["law"]=="mono":
        phi0 = np.array([phi["angle"]] * Nconf)
      elif phi["law"]=="iso":
        try:
          mangle = theta["max"]*np.pi/180.
        except KeyError:
          mangle = np.pi
        phi0 = np.random.uniform(0,mangle,Nconf)
      elif phi["law"]=="gauss":
        phi0 = np.random.normal(phi["mu"],phi["sigma"],Nconf)
      # Generate energy
      if ekin["law"]=="mono":
        ekin0 = np.array([ekin["E"]] * Nconf)
      elif ekin["law"]=="exp":
        ekin0 = np.random.exponential(ekin["T"],Nconf)

      # Reconstruct momentum from energy and angle distributions
      mass = self._ps.mass
      etot = ekin0 + mass
      px = np.sqrt((etot**2 - mass**2)/(1. + np.tan(theta0)**2 + np.tan(phi0)**2))
      py = px * np.tan(theta0)
      pz = px * np.tan(phi0)

      # Generate position
      if pos is None:
        x = np.array([0.] * Nconf)
        y = np.array([0.] * Nconf)
        z = np.array([0.] * Nconf)

      # Generate time
      if time is None:
        t = np.array([0.] * Nconf)

      if verbose: print("Done !")
      # Update current object
      self.update(w,x,y,z,px,py,pz,t,verbose=verbose)

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

  def lorrentz(self,beta_CM,pos_CM=None):
    """
    Lorrentz-transformate the particle phase-space with given speed of the center of mass.
    """
    # If no position given, default is the origin
    if pos_CM is None:
      pos_CM = [0.,0.,0.]
    # Calculate gamma factor from beta_CM
    # beta_CM = np.array(beta_CM)
    # v_CM = beta_CM * 2.99792458e8 * 1e6/1e15 # speed of CM in um/fs
    # gamma_CM = 1./np.sqrt(1-beta_CM**2)
    # # Time transformation
    # t = gamma_CM*(self.t
    # # Position transformation
    # # x = (self.x - pos_CM[0])*gamma_CM[0]
    # # y = (self.y - pos_CM[1])*gamma_CM[1]
    # # z = (self.z - pos_CM[2])*gamma_CM[2]
    # x = self.x
    # y = self.y
    # z = self.z
    # # Momentum transformation
    # px = self.px*gamma_CM[0]
    # py = self.px*gamma_CM[1]
    # pz = self.px*gamma_CM[2]
    # # No transformations on statistical weights
    # w = self.w

    # beta_CM = np.array(beta_CM)
    # b = np.dot(beta_CM,beta_CM)
    # n = beta_CM/b
    # c = 2.99792458e8 * 1e6/1e15 # speed of CM in um/fs
    # v = b * c
    # gamma = 1./np.sqrt(1-beta_CM**2)
    # g = 1./np.sqrt(1-b**2)
    # G = np.matrix([
    # [g           ,-g*b*n[0]          ,-g*b*n[1]          ,-g*b*n[2]      ],
    # [-g*b*n[0]   ,1+(g-1)*n[0]**2    ,(g-1)*n[0]*n[1]    ,(g-1)*n[0]*n[2]],
    # [-g*b*n[1]   ,(g-1)*n[1]*n[0]    ,1+(g-1)*n[1]**2    ,(g-1)*n[1]*n[2]],
    # [-g*b*n[2]   ,(g-1)*n[2]*n[0]    ,(g-1)*n[2]*n[1]    ,1+(g-1)*n[2]**2]
    # ])
    # print(G)
    # # Time transformation
    # q1 = [c*self.t.copy(),self.x.copy(),self.y.copy(),self.z.copy()]
    # p1 = [self.ekin.copy()/c,self.px.copy(),self.py.copy(),self.pz.copy()]
    # q2 = np.dot(G,q1)
    # p2 = np.dot(G,p1)
    # w = self.w
    # t = np.array(q2[0]/c)[0]
    # x = np.array(q2[1])[0]
    # y = np.array(q2[2])[0]
    # z = np.array(q2[3])[0]
    # px = np.array(p1[1])[0]
    # py = np.array(p1[2])[0]
    # pz = np.array(p1[3])[0]
    # self.update(w,x,y,z,px,py,pz,t)

    beta_CM = np.array(beta_CM)
    b = np.dot(beta_CM,beta_CM)
    n = beta_CM/b
    c = 2.99792458e8 * 1e6/1e15 # speed of CM in um/fs
    v = beta_CM * c
    gamma = 1./np.sqrt(1-beta_CM**2)
    g = 1./np.sqrt(1-b**2)

    w1 = self.w
    t1 = self.t
    r1 = np.array([self.x,self.y,self.z])
    p1 = np.array([self.px,self.py,self.pz])

    w2 = w1
    t2 = g*(t1 + sum(r1*v)/c**2)
    r2 = r1 + ((g-1)/np.dot(v,v) * sum(r1*v) + g*t1)*v
    p2 = p1

    self.update(w2,r2[0],r2[1],r2[2],p2[0],p2[1],p2[2],t2)


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
