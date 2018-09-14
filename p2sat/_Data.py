#coding:utf8
import numpy as np

class _Data(object):
  """
  Class containing raw data and methods to manipulate it.

  Attributes
  ----------
  w : numpy.ndarray
    statistical weight
  x,y,z : numpy.ndarray
    x,y,z position in um
  px,py,pz : numpy.ndarray
    momentum in x,y,z direction in MeV/c
  t : numpy.ndarray
    time in fs
  r : numpy.ndarray
    absolute distance to the x axis in um
  p : numpy.ndarray
    absolute momentum in MeV/c
  ekin : numpy.ndarray
    energy in MeV
  gamma : numpy.ndarray
    gamma factor
  beta : numpy.ndarray
    normalized velocity
  v : numpy.ndarray
    velocity in um/fs
  theta : numpy.ndarray
    polar angle (between px and the plane (py,pz)) in degree
  phi : numpy.ndarray
    azimutal angle (between py and pz) in degree

  Notes
  -----
  As all the calculations are done with the previously defined units,
  the input data might be firstly converted to those units.

  Calculations :

  - r is defined as :math:`\sqrt{y^2+z^2}`
  - p is defined as :math:`\sqrt{p_x^2+p_y^2+p_z^2}`
  - theta is defined as :math:`\\arccos{p_x/p}`
  - phi is defined as :math:`\\arctan{p_z/p_y}`
  - ekin is defined as

    - :math:`(\sqrt{(p/m_e c)^2+1}-1) \\times m_e c^2` for massive species
    - :math:`p` otherwise (here ekin is the total particle energy)
  - gamma is defined as

    - :math:`E_{kin}/m_e c^2 + 1` for massive species
    - :math:`+\infty` otherwise

  References
  ----------
  For ekin, gamma, beta, v :
  https://en.wikipedia.org/wiki/Lorentz_factor

  For theta, phi :
  https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
  with switching (x->py, y->pz, z->px)
  A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system
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
              'beta'  : '$\\beta$',
              'v'     : '$v$',

              'theta' : '$\\theta$',
              'phi'   : '$\\phi$',

              'e-'    : '$e^-$',
              'e+'    : '$e^+$',
              'mu-'   : '$\mu^-$',
              'mu+'   : '$\mu^+$'
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
              'beta'  : None,
              'v'     : 'um/fs',

              'theta' : 'deg',
              'phi'   : 'deg'
              }

  def update(self,w,x,y,z,px,py,pz,t,verbose=True):
    """
    Update class attributes with new values.

    Parameters
    ----------
    w,x,y,z,px,py,pz,t : list or numpy.ndarray
      particle phase space. More information can be found in data object documentation
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
    self.theta  = np.degrees(np.arccos(self.px/self.p))
    self.phi    = np.degrees(np.arctan2(self.pz,self.py))
    mass = self._ps.mass
    c = 2.99792458e8 * 1e6/1e15 # speed of light in um/fs
    if mass == 0:
      self.ekin   = self.p
      self.gamma  = np.array([np.inf]*len(w))
      self.beta   = np.array([1.]*len(w))
      self.v      = np.array([c]*len(w))
    else:
      self.ekin   = (np.sqrt((self.p/mass)**2 + 1) - 1) * mass
      self.gamma  = self.ekin/mass + 1.
      self.beta   = np.sqrt(1.-1/self.gamma**2)
      self.v      = self.beta * c
    if verbose: print("Done !")

  def generate(self,Nconf,Npart,ekin,theta,phi,x=None,y=None,z=None,r=None,t=None,verbose=True):
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

      Details of the calculations :

      Considering :math:`E_T = E_k + E_m` being the total energy, with
      :math:`E_k` the kinetic energy and :math:`E_m` the rest mass energy, we have
      :math:`E_T^2 = p^2 + E_m^2` and :math:`p^2=p_x^2 + p_y^2 + p_z^2` so
      :math:`p_x^2+p_y^2+p_z^2=E_T^2 - E_m^2`.

      Assuming :math:`\\tan{\\theta}=\\frac{p_y}{p_x}` and
      :math:`\\tan{\\phi}=\\frac{p_z}{p_x}` we get

      :math:`p_x = \sqrt{\\frac{E_T^2 - E_m^2}{1 + \\tan{\\theta}^2 + \\tan{\phi}^2}}`
      :math:`p_y = p_x \\tan{\\theta}`
      :math:`p_z = p_x \\tan{\phi}`

      Examples
      --------
      Assuming a `PhaseSpace` object is instanciated as `eps`, you can generate
      a mono-energetic source in isotropic direction for 1e12 particles represented
      by 1e6 configurations as follows

      >>> eps.data.generate(Nconf=1e6,Npart=1e12,ekin={"law":"mono","E":20.0},theta={"law":"iso"},phi={"law":"iso"})
      """
      # Print a starting message
      if verbose:
        print("Generate %s phase-space for \"%s\" ekin law, \"%s\" theta law, \"%s\" phi law ..."
        %(self._ps.specie,ekin["law"],theta["law"],phi["law"]))

      # Ensure that Nconf is of type int (for values such as 1e6)
      Nconf = int(Nconf)
      # Generate weights
      weight = float(Npart)/Nconf
      w = np.array([weight] * Nconf)

      # Generate theta angle
      if theta["law"]=="mono":
        theta0 = np.array([theta["theta0"]] * Nconf)
      elif theta["law"]=="iso":
        try:
          mangle = np.radians(theta["max"])
        except KeyError:
          mangle = np.pi
        theta0 = np.random.uniform(0.,mangle,Nconf)
      elif theta["law"]=="gauss":
        mu = np.radians(theta["mu"])
        sigma = np.radians(theta["sigma"])
        theta0 = abs(np.random.normal(mu,sigma,Nconf))
      # Generate phi angle
      if phi["law"]=="mono":
        phi0 = np.array([phi["theta0"]] * Nconf)
      elif phi["law"]=="iso":
        try:
          mangle = np.radians(phi["max"])
        except KeyError:
          mangle = np.pi
        phi0 = np.random.uniform(-mangle,mangle,Nconf)
      elif phi["law"]=="gauss":
        mu = np.radians(phi["mu"])
        sigma = np.radians(phi["sigma"])
        phi0 = np.random.normal(mu,sigma,Nconf)
      # Generate energy
      if ekin["law"]=="mono":
        ekin0 = np.array([ekin["ekin0"]] * Nconf)
      elif ekin["law"]=="exp":
        ekin0 = np.random.exponential(ekin["ekin0"],Nconf)

      # Reconstruct momentum from energy and angle distributions
      mass  = self._ps.mass
      p     = np.sqrt(ekin0**2 + 2*ekin0*mass)
      px    = p * np.cos(theta0)
      py    = np.sign(phi0)*np.sqrt((p**2 - px**2)/(1. + np.tan(phi0)**2))
      pz    = py*np.tan(phi0)

      # Generate position
      if x is None:
        x0 = np.array([0.] * Nconf)
      elif x["law"]=="mono":
        x0 = np.array([x["x0"]] * Nconf)
      elif x["law"]=="exp":
        x0 = np.random.exponential(x["x0"],Nconf)
      elif x["law"]=="gauss":
        mu = np.radians(x["mu"])
        sigma = np.radians(x["sigma"])
        x0 = np.random.normal(mu,sigma,Nconf)

      if y is None:
        y0 = np.array([0.] * Nconf)
      elif y["law"]=="mono":
        y0 = np.array([y["y0"]] * Nconf)
      elif y["law"]=="exp":
        y0 = np.random.exponential(y["y0"],Nconf)
      elif y["law"]=="gauss":
        mu = np.radians(y["mu"])
        sigma = np.radians(y["sigma"])
        y0 = np.random.normal(mu,sigma,Nconf)

      if z is None:
        z0 = np.array([0.] * Nconf)
      elif z["law"]=="mono":
        z0 = np.array([z["z0"]] * Nconf)
      elif z["law"]=="exp":
        z0 = np.random.exponential(z["z0"],Nconf)
      elif z["law"]=="gauss":
        mu = np.radians(z["mu"])
        sigma = np.radians(z["sigma"])
        z0 = np.random.normal(mu,sigma,Nconf)

      if r is None:
        r0 = np.array([0.] * Nconf)
      elif r["law"]=="gauss":
        mu = r["mu"]
        sigma = r["sigma"]
        r0 = np.random.normal(mu,sigma,Nconf)
        angle = np.random.uniform(0.,2*np.pi,Nconf)
        z0 = np.sign(r0) * np.sqrt(r0**2/(1. + np.tan(angle)**2))
        y0 = z0 * np.tan(angle)

      # Generate time
      if t is None:
        t0 = np.array([0.] * Nconf)
      elif t["law"]=="mono":
        t0 = np.array([t["t0"]] * Nconf)
      elif t["law"]=="exp":
        t0 = np.random.exponential(t["t0"],Nconf)
      elif t["law"]=="gauss":
        mu = np.radians(t["mu"])
        sigma = np.radians(t["sigma"])
        t0 = np.random.normal(mu,sigma,Nconf)

      if verbose: print("Done !")
      # Update current object
      self.update(w,x0,y0,z0,px,py,pz,t0,verbose=verbose)

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

  def propagate(self,x_pos=None,time=None,verbose=True):
    """
    Propagate the phase space to a given position or time.

    Parameters
    ----------
    x_pos : float, optional
      propagate the phase-space untill x = x_pos. Default is None (no propagation)
    time : float, optional
      propagate the phase-space untill t = time. Default is None (no propagation)

    Notes
    -----
    x_pos and time can not be defined simultaneously.
    """
    if time is None and x_pos is None:
      raise ValueError("You must define time or x_pos.")
    if time is not None and x_pos is not None:
      raise ValueError("time and x_pos can not be defined simultaneously.")

    w = self.w
    px = self.px
    py = self.py
    pz = self.pz

    if time is not None:
      if verbose: print("Propagate %s phase-space to time = %.4E fs."%(self._ps.specie,time))
      t = np.array([time]*len(w))
      Dt = t - self.t
      x = self.x + (self.px/self.p)*self.v*Dt
      y = self.y + (self.py/self.p)*self.v*Dt
      z = self.z + (self.pz/self.p)*self.v*Dt

    if x_pos is not None:
      if verbose: print("Propagate %s phase-space to x = %.4E um."%(self._ps.specie,x_pos))
      x = np.array([x_pos]*len(w))
      Dt = (x - self.x)/self.v
      t = self.t + Dt
      y = self.y + (self.py/self.p)*self.v*Dt
      z = self.z + (self.pz/self.p)*self.v*Dt

    self.update(w,x,y,z,px,py,pz,t,verbose=verbose)

  def lorentz(self,beta_CM,verbose=True):
    """
    Lorentz-transformate the particle phase-space with given speed of the center of mass.

    Notes
    -----
    https://en.wikipedia.org/wiki/Lorentz_transformation#Transformation_of_other_quantities
    """
    if verbose: print("Lorentz-transform %s phase-space with center of mass frame moving at %s c."%(self._ps.specie,beta_CM))
    # lowercase : scalar, caption : vector

    B = -np.array(beta_CM)
    b = np.dot(B,B)
    N = B/b
    c = 2.99792458e8 * 1e6/1e15 # speed of light in um/fs
    v = b * c
    g = 1./np.sqrt(1-b**2)

    # Four position
    a1 = c*self.t
    Z1 = np.array([self.x,self.y,self.z]).T

    a2,Z2 =[],[]
    for i in range(len(a1)):
        a2.append(g*(a1[i] - (v*np.dot(N,Z1[i]))/c))
        Z2.append(Z1[i] + (g-1)*np.dot(Z1[i],N)*N - g*(a1[i]*v*N)/c)

    t = np.array(a2)/c
    x,y,z = np.array(Z2).T

    # Four momentum
    a1 = self.ekin/c
    Z1 = np.array([self.px,self.py,self.pz]).T

    a2,Z2 =[],[]
    for i in range(len(a1)):
        a2.append(g*(a1[i] - (v*np.dot(N,Z1[i]))/c))
        Z2.append(Z1[i] + (g-1)*np.dot(Z1[i],N)*N - g*(a1[i]*v*N)/c)

    px,py,pz = np.array(Z2).T

    w = self.w
    self.update(w,x,y,z,px,py,pz,t)

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

    if with_time:
      bi,hi=hn([self.x,self.y,self.z,self.px,self.py,self.pz,self.t],normed=False,**kargs)
      bx,by,bz,bpx,bpy,bpz,bt=bi
    else:
      bi,hi=hn([self.x,self.y,self.z,self.px,self.py,self.pz],normed=False,**kargs)
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
