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
        kinetic energy (for massive species) or total energy (for massless species) in MeV
    gamma : numpy.ndarray
        Lorentz factor
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

    Definitions :

    - :math:`r = \sqrt{y^2+z^2}`
    - :math:`p = \sqrt{p_x^2+p_y^2+p_z^2}`
    - :math:`\\theta = \\arccos{p_x/p}`
    - :math:`\phi = \\arctan{p_z/p_y}`

    - For massive species (with mass :math:`m`)

        - :math:`E_{kin} = (\sqrt{(p/m c)^2+1}-1) \\times m c^2`
        - :math:`\gamma = E_{kin}/m c^2 + 1`
        - :math:`\\beta = \sqrt{1 - \\frac{1}{\gamma^2}}`
        - :math:`v = \\beta \\times c`

    - For massless species

        - :math:`E_{kin} = p`
        - :math:`\gamma = +\infty`
        - :math:`\\beta = 1`
        - :math:`v = c`

    All the attributes can not be overwriden as they are defined as properties.
    Please call the `update` method to update particle phase-space data.

    References
    ----------
    For ekin, gamma, beta, v :
    https://en.wikipedia.org/wiki/Lorentz_factor

    For theta, phi :
    https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
    with switching (x->py, y->pz, z->px).
    A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system
    """
    def __init__(self,PhaseSpace):
        self._ps= PhaseSpace
        self.update(None,None,None,None,None,None,None,None,verbose=False)

        specie_label = self._ps.specie['label']
        if self._ps.specie['name']=="gamma":
            etot_label = 'E_{\gamma}'
            ekin_label = 'E_{\gamma}'
        else:
            etot_label = 'E_{Tot}'
            ekin_label = 'E_{kin}'

        self.labels = {}                ; self.units = {}
        l = self.labels                 ; u = self.units

        l['w'] = 'N_{%s}'%specie_label  ; u['w'] = None
        l['x'] = 'x'                    ; u['x'] = u'µm'
        l['y'] = 'y'                    ; u['y'] = u'µm'
        l['z'] = 'z'                    ; u['z'] = u'µm'
        l['px'] = 'p_x'                 ; u['px'] = 'MeV/c'
        l['py'] = 'p_y'                 ; u['py'] = 'MeV/c'
        l['pz'] = 'p_z'                 ; u['pz'] = 'MeV/c'
        l['t'] = 't'                    ; u['t'] = 'fs'

        l['r'] = 'r'                    ; u['r'] = u'um'
        l['p'] = 'p'                    ; u['p'] = 'MeV/c'

        l['etot'] = etot_label          ; u['etot'] = 'MeV'
        l['ekin'] = ekin_label          ; u['ekin'] = 'MeV'
        l['beta'] = '\\beta'            ; u['beta'] = None
        l['gamma'] = '\gamma'           ; u['gamma'] = None
        l['v'] = 'v'                    ; u['v'] = u'um/fs'

        l['theta'] = '\\theta'          ; u['theta'] = 'deg'
        l['phi'] = '\\phi'              ; u['phi'] = 'deg'
        l['omega'] = '\Omega'           ; u['omega'] = 'sr'

        l['ekin_density'] = '(N_{%s} %s)'%(specie_label,ekin_label)
        u['ekin_density'] = 'MeV'

    @property
    def w(self):
        return self._w

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @property
    def px(self):
        return self._px

    @property
    def py(self):
        return self._py

    @property
    def pz(self):
        return self._pz

    @property
    def t(self):
        return self._t

    @property
    def r(self):
        return np.sqrt(self.y**2+self.z**2)

    @property
    def p(self):
        return np.sqrt(self.px**2+self.py**2+self.pz**2)

    @property
    def theta(self):
        return np.degrees(np.arccos(self.px/self.p))

    @property
    def phi(self):
        return np.degrees(np.arctan2(self.pz,self.py))

    @property
    def omega(self):
        return np.pi*self.r**2/self.x**2

    @property
    def etot(self):
        mass = self._ps.specie["mass"]
        return np.sqrt(self.p**2 + mass**2)

    @property
    def ekin(self):
        mass = self._ps.specie["mass"]
        return self.etot-mass

    @property
    def gamma(self):
        mass = self._ps.specie["mass"]
        return self.ekin/mass + 1.

    @property
    def beta(self):
        return np.sqrt(1.-1./self.gamma**2)

    @property
    def v(self):
        c = 2.99792458e8 * 1e6/1e15 # speed of light in um/fs
        return self.beta * c

    @property
    def ekin_density(self):
        return self.ekin * self.w

    # @property
    # def brilliance(self):
    #     return photons/s/mm2/mrad2/0.1%BW.

    def get_ps(self):
        """
        Return current phase space raw data.
        """
        return (self.w,self.x,self.y,self.z,self.px,self.py,self.pz,self.t)

    def get_axis(self,axis,select=None):
        """
        Return given axis

        Parameters
        ----------
        axis : str
            axis name
        select : dict
            filtering dictionnary

        Examples
        --------

        >>> eps.data.get_axis("w")
        >>> eps.data.get_axis("w",select={'x':150,'ekin':[0.511,None]})
        """
        if type(axis) is not str:
            # If not a str, return axis
            return axis
        else:
            # Evaluate axis
            ax = eval("self.%s"%axis)
            # Filter the data if needed
            if select is not None:
                ax = self.select(ax,faxes=select.keys(),frange=select.values())
            # return result
            return ax


    def update(self,w,x,y,z,px,py,pz,t,verbose=True):
        """
        Update class attributes with new values.

        Parameters
        ----------
        w,x,y,z,px,py,pz,t : list or numpy.ndarray
          particle phase space. More information can be found in data object documentation
        verbose : bool
          verbosity of the function. If True, a message is displayed when the attributes are loaded in memory
        """
        if verbose: print("Updating raw values ...")
        # Save values into np.array objects
        self._w  = np.array(w)
        self._x  = np.array(x)
        self._y  = np.array(y)
        self._z  = np.array(z)
        self._px = np.array(px)
        self._py = np.array(py)
        self._pz = np.array(pz)
        self._t  = np.array(t)
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
        x,y,z : dict, optional
            parameters to generate position distribution. Default is 0 for x,y,z
        r : dict,optional
            parameters to generate transverse position distribution. Default is 0
        t : dict, optional
            parameters to generate time distribution. Default is 0

        Notes
        -----
        The dictionnaries must each time at least contain the key 'law' with a value
        depending on which law are available for each physical quantity

        For dict `ekin`, available laws are :

        - 'mono', for a mono-energetic source. Energy must be given as a value of keyword 'ekin0'
        - 'exp', for exponential energy. Charasteristic energy must be given as a value of keyword 'ekin0'

        For dict `theta` and `phi`, available laws are :

        - 'mono', for a directional source. Angle must be given as a value of keyword 'theta0' or '/phi0'
        - 'iso', for an isotropic source. Optional keywords 'min' and 'max' can be given to specify a minimum/maximum angle (centered on 0 deg)
        - 'gauss', for a gaussian spreading. Center of the distribution must be given with keyword 'mu', and standard deviantion with keyword 'sigma'

        For dict `x`,`y`,`z`,`t`, available laws are :

        - 'mono', for a unique position/time. This parameter must be given as a value of keyword `x0`/`y0`/`z0`/`t0`
        - 'range', for a uniform law in a given range. Keywords 'min' and 'max' MUST be given to specify a minimum and maximum
        - 'exp', for exponential distribution. Charasteristic length/time must be given as a value of keyword `x0`/`y0`/`z0`/`t0`
        - 'gauss', for a gaussian distribution. Center of the distribution must be given with keyword 'mu', and standard deviantion with keyword 'sigma'

        For dict `r`, available laws are :

        - 'range', for a uniform law in a given range. Keywords 'min' and 'max' MUST be given to specify a minimum and maximum
        - 'gauss', for a gaussian distribution. Center of the distribution must be given with keyword 'mu', and standard deviantion with keyword 'sigma'

        Details of the calculations :

        Considering :math:`E_T = E_k + E_m` being the total energy, with
        :math:`E_k` the kinetic energy and :math:`E_m` the rest mass energy.

        We also have :math:`E_T^2 = p^2 + E_m^2` and :math:`p^2=p_x^2 + p_y^2 + p_z^2`
        with :math:`p` in MeV/c.

        Assuming :math:`\cos{\\theta}=\\frac{p_x}{p}` and
        :math:`\\tan{\\phi}=\\frac{p_z}{p_y}` we finaly get

        - :math:`p = E_k^2 - 2 E_k E_m`
        - :math:`p_x = p \cos{\\theta}`
        - :math:`p_y = \\frac{\phi}{\mid \phi \mid} \sqrt{\\frac{p^2 - p_x^2}{1 + \\tan^2{\phi}}}`
        - :math:`p_z = p_y \\tan{\phi}`

        Examples
        --------
        With a `PhaseSpace` object instanciated as `eps`, you can generate
        a mono-energetic source in isotropic direction for 1e12 particles represented
        by 1e6 configurations as follows

        >>> eps.data.generate(Nconf=1e6,Npart=1e12,
        ...                  ekin={"law":"mono","ekin0":20.0},
        ...                  theta={"law":"iso"},
        ...                  phi={"law":"iso"})
        ...
        """
        # Print a starting message
        if verbose:
            print("Generate %s phase-space ..."%(self._ps.specie["name"]))
            print("    ekin  : %s"%ekin)
            print("    theta : %s"%theta)
            print("    phi   : %s"%phi)
            if x is not None: print("    x     : %s"%x)
            if y is not None: print("    y     : %s"%y)
            if z is not None: print("    z     : %s"%z)
            if r is not None: print("    r     : %s"%r)
            if t is not None: print("    t     : %s"%t)

        # Ensure that Nconf is of type int (for values such as 1e6)
        Nconf = int(Nconf)
        # Generate weights
        weight = float(Npart)/Nconf
        g_w = np.array([weight] * Nconf)

        # Generate theta angle
        if theta["law"]=="mono":
            g_theta = np.array([theta["theta0"]] * Nconf)
        elif theta["law"]=="iso":
            try:
                mini = np.radians(theta["min"])
            except KeyError:
                mini = 0.
            try:
                maxi = np.radians(theta["max"])
            except KeyError:
                maxi = np.pi
            g_theta = np.random.uniform(0.,maxi,Nconf)
        elif theta["law"]=="gauss":
            mu = np.radians(theta["mu"])
            sigma = np.radians(theta["sigma"])
            g_theta = abs(np.random.normal(mu,sigma,Nconf))
        # Generate phi angle
        if phi["law"]=="mono":
            g_phi = np.array([phi["phi0"]] * Nconf)
        elif phi["law"]=="iso":
            try:
                mini = np.radians(phi["min"])
            except KeyError:
                mini = 0.
            try:
                maxi = np.radians(phi["max"])
            except KeyError:
                maxi = np.pi
            g_phi = np.random.uniform(-maxi,maxi,Nconf)
        elif phi["law"]=="gauss":
            mu = np.radians(phi["mu"])
            sigma = np.radians(phi["sigma"])
            g_phi = np.random.normal(mu,sigma,Nconf)
        # Generate energy
        if ekin["law"]=="mono":
            g_ekin = np.array([ekin["ekin0"]] * Nconf)
        elif ekin["law"]=="exp":
            g_ekin = np.random.exponential(ekin["ekin0"],Nconf)

        # Reconstruct momentum from energy and angle distributions
        mass  = self._ps.specie["mass"]
        g_p     = np.sqrt(g_ekin**2 + 2*g_ekin*mass)
        g_px    = g_p * np.cos(g_theta)
        g_py    = np.sign(g_phi)*np.sqrt((g_p**2 - g_px**2)/(1. + np.tan(g_phi)**2))
        g_pz    = g_py*np.tan(g_phi)

        # Generate position
        if x is None:
            g_x = np.array([0.] * Nconf)
        elif x["law"]=="mono":
            g_x = np.array([x["x0"]] * Nconf)
        elif x["law"]=="range":
            g_x = np.random.uniform(x["min"],x["max"],Nconf)
        elif x["law"]=="exp":
            g_x = np.random.exponential(x["x0"],Nconf)
        elif x["law"]=="gauss":
            mu = x["mu"]
            sigma = x["sigma"]
            g_x = np.random.normal(mu,sigma,Nconf)

        if y is None:
            g_y = np.array([0.] * Nconf)
        elif y["law"]=="mono":
            g_y = np.array([y["y0"]] * Nconf)
        elif y["law"]=="range":
            g_y = np.random.uniform(y["min"],y["max"],Nconf)
        elif y["law"]=="exp":
            g_y = np.random.exponential(y["y0"],Nconf)
        elif y["law"]=="gauss":
            mu = y["mu"]
            sigma = y["sigma"]
            g_y = np.random.normal(mu,sigma,Nconf)

        if z is None:
            g_z = np.array([0.] * Nconf)
        elif z["law"]=="mono":
            g_z = np.array([z["z0"]] * Nconf)
        elif z["law"]=="range":
            g_z = np.random.uniform(z["min"],z["max"],Nconf)
        elif z["law"]=="exp":
            g_z = np.random.exponential(z["z0"],Nconf)
        elif z["law"]=="gauss":
            mu = z["mu"]
            sigma = z["sigma"]
            g_z = np.random.normal(mu,sigma,Nconf)

        if r is None:
            g_r = np.array([0.] * Nconf)
        elif r["law"]=="range":
            g_r = np.random.uniform(r["min"],r["max"],Nconf)
            angle = np.random.uniform(0.,2*np.pi,Nconf)
            g_z = g_r * np.cos(angle)
            g_y = g_r * np.sin(angle)
        elif r["law"]=="gauss":
            mu = r["mu"]
            sigma = r["sigma"]
            g_r = np.random.normal(mu,sigma,Nconf)
            angle = np.random.uniform(0.,2*np.pi,Nconf)
            g_z = g_r * np.cos(angle)
            g_y = g_r * np.sin(angle)

        # Generate time
        if t is None:
            g_t = np.array([0.] * Nconf)
        elif t["law"]=="mono":
            g_t = np.array([t["t0"]] * Nconf)
        elif t["law"]=="range":
            g_t = np.random.uniform(t["min"],t["max"],Nconf)
        elif t["law"]=="exp":
            g_t = np.random.exponential(t["t0"],Nconf)
        elif t["law"]=="gauss":
            mu = t["mu"]
            sigma = t["sigma"]
            g_t = np.random.normal(mu,sigma,Nconf)

        if verbose: print("Done !")
        # Update current object
        self.update(g_w,g_x,g_y,g_z,g_px,g_py,g_pz,g_t,verbose=verbose)

    def select(self,axis,faxes,frange,fpp=1e-7):
        """
        Filter an axis with a value/range on another axis.

        Parameters
        ----------
        axis : str or numpy.ndarray
            axis to filter
        faxes : list of str or list of numpy.ndarray
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
        Given the `PhaseSpace` instance `eps`, it is possible to filter by an int value
        (Select all the :math:`w` satisfying :math:`x=3`)

        >>> w,x = eps.data.select('w','x',3)

        or filter by a range (Select all the :math:`\\theta` with :math:`E_{kin} \in [0.511,+\infty]` MeV)

        >>> theta,ekin = eps.data.select('theta','ekin',[0.511,None])

        If frange is a list/tuple or a float, the filtering is done with a fpp precision
        """
        # Get a copy of axis and faxes (from a str or not)
        axis=self.get_axis(axis,select=None)
        for i,fax in enumerate(faxes):
            faxes[i] = self.get_axis(fax,select=None)

        # Loop over filtering axes
        for i,_ in enumerate(faxes):
            # Filter in a given range
            if isinstance(frange[i],(list,tuple,type(np.array(0)))):
                # Default ranges are min and max of filtering axis
                if frange[i][0] is None: frange[i][0]=min(faxes[i])
                if frange[i][1] is None: frange[i][1]=max(faxes[i])
                # The filtr array is an array of bool, whether the filtering condition is fullfiled or not
                filtr=np.array([x>frange[i][0]*(1-fpp) and x<frange[i][1]*(1+fpp) for x in faxes[i]])
                axis=axis[filtr]
            # Filter for an int
            elif type(frange[i]) is int:
                axis=axis[faxes[i]==frange[i]]
            # Filter for a float
            elif type(frange[i]) is float:
                axis=self.select(axis,faxes=[faxes[i]],frange=[[frange[i],frange[i]]])
            else:
                raise TypeError('frange type must be int/float or list/tuple of 2 float.')

            # filter next faxes with current faxes
            if len(faxes)>i+1:faxes[i+1]=self.select(faxes[i+1],[faxes[i]],[frange[i]])

        return axis

    def full_select(self,faxes,frange,fpp=1e-7,update=False,verbose=True):
        """
        Select all the phase space with given condition

        Parameters
        ----------
        faxes : list of str or list of numpy.ndarray
            filtering axis
        frange : list of int, float, list/tuple of 2 float
            filtering value/range (value if int, range if float or list/tuple). If a frange element is None, the minimum/maximum value is taken
        fpp : float, optional
            relative floating point precision. Default is 1e-7
        update : bool, optional
            update or not the current `PhaseSpace` instance. Default is false
        verbose : bool, optional
            verbosity

        See Also
        --------
        data.select
        """
        if verbose: print("Filtering %s phase space with axis %s ..."%(self._ps.specie["name"],faxes))
        data = []
        for ax in self.get_ps():
            data.append(self.select(ax,faxes,frange,fpp=fpp))

        w   = data[0::8]
        x   = data[1::8]
        y   = data[2::8]
        z   = data[3::8]
        px  = data[4::8]
        py  = data[5::8]
        pz  = data[6::8]
        t   = data[7::8]

        if verbose: print("Done !")

        if update:
            self.update(w,x,y,z,px,py,pz,t,verbose=verbose)
        else:
            return w,x,y,z,px,py,pz,t

    def transformate(self,T=None,R=None,rotate_first=False,verbose=True):
        """
        Transformate the particle phase space with given translation and rotation.

        Parameters
        ----------
        T : tuple of 3 float, optional
            translate (x,y,z) position of T um. Default is (0,0,0)
        R : tuple of 3 float, optional
            rotate (x,y,z) and (px,py,pz) of R degree. Default is (0,0,0)
        rotate_first : bool, optional
            process rotation before translation. Default is false
        verbose : bool, optional
            verbosity
        """
        # Set defaults
        if T is None: T = (0.,0.,0.)
        if R is None: R = (0.,0.,0.)
        if verbose:
            print("Transformate %s phase space ..."%self._ps.specie["name"])
            print("    translation : (%.2E,%.2E,%.2E)"%(T[0],T[1],T[2]))
            print("    rotation    : (%.2E,%.2E,%.2E)"%(R[0],R[1],R[2]))

        # Define translation matrix
        Tx = T[0]
        Ty = T[1]
        Tz = T[2]

        # Define rotation matrix
        a = np.radians(R[0])
        Rx = np.array([
            [1              ,0              ,0              ],
            [0              ,np.cos(a)      ,-np.sin(a)     ],
            [0              ,np.sin(a)      ,np.cos(a)      ]])
        a = np.radians(R[1])
        Ry = np.array([
            [np.cos(a)      ,0              ,np.sin(a)      ],
            [0              ,1              ,0              ],
            [-np.sin(a)     ,0              ,np.cos(a)      ]])
        a = np.radians(R[2])
        Rz = np.array([
            [np.cos(a)      ,-np.sin(a)     ,0              ],
            [np.sin(a)      ,np.cos(a)      ,0              ],
            [0              ,0              ,1              ]])

        # Get current phase space data
        w,x,y,z,px,py,pz,t = self.get_ps()

        # Translate position
        if not rotate_first:
            x += Tx
            y += Ty
            z += Tz

        X,Y,Z,Px,Py,Pz=[],[],[],[],[],[]
        R = np.dot(np.dot(Rx,Ry),Rz)
        # Loop over all configurations
        for i,_ in enumerate(w):
            # Rotate position
            rX = np.dot([x[i],y[i],z[i]],R)
            X.append(rX[0])
            Y.append(rX[1])
            Z.append(rX[2])
            # Rotate momentum
            rP = np.dot([px[i],py[i],pz[i]],R)
            Px.append(rP[0])
            Py.append(rP[1])
            Pz.append(rP[2])

        if rotate_first:
            X = np.array(X) + Tx
            Y = np.array(Y) + Ty
            Z = np.array(Z) + Tz

        # Update raw data
        self.update(w,X,Y,Z,Px,Py,Pz,t,verbose=verbose)

    def propagate(self,x_pos=None,time=None,update=True,verbose=True):
        """
        Propagate the phase space to a given position or time.

        Parameters
        ----------
        x_pos : float, optional
            propagate the phase-space untill x = x_pos. Default is None (no propagation)
        time : float, optional
            propagate the phase-space untill t = time. Default is None (no propagation)
        verbose : bool, optional
            verbosity

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
            if verbose: print("Propagate %s phase-space to time = %.4E fs."%(self._ps.specie["name"],time))
            t = np.array([time]*len(w))
            Dt = t - self.t
            x = self.x + (self.px/self.p)*self.v*Dt
            y = self.y + (self.py/self.p)*self.v*Dt
            z = self.z + (self.pz/self.p)*self.v*Dt

        if x_pos is not None:
            if verbose: print("Propagate %s phase-space to x = %.4E um."%(self._ps.specie["name"],x_pos))
            x = np.array([x_pos]*len(w))
            Dt = (x - self.x)/self.v
            t = self.t + Dt
            y = self.y + (self.py/self.p)*self.v*Dt
            z = self.z + (self.pz/self.p)*self.v*Dt

        if update:
            self.update(w,x,y,z,px,py,pz,t,verbose=verbose)
        else:
            return w,x,y,z,px,py,pz,t

    def lorentz(self,beta_CM,verbose=True):
        """
        Lorentz-transformate the particle phase-space with given speed of the center of mass.

        TODO

        References
        ----------
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

        # 4 position
        a1 = c*self.t
        Z1 = np.array([self.x,self.y,self.z]).T

        a2,Z2 =[],[]
        for i in range(len(a1)):
            a2.append(g*(a1[i] - (v*np.dot(N,Z1[i]))/c))
            Z2.append(Z1[i] + (g-1)*np.dot(Z1[i],N)*N - g*(a1[i]*v*N)/c)

        t = np.array(a2)/c
        x,y,z = np.array(Z2).T

        # 4 momentum
        a1 = self.ekin/c
        Z1 = np.array([self.px,self.py,self.pz]).T

        a2,Z2 =[],[]
        for i in range(len(a1)):
            a2.append(g*(a1[i] - (v*np.dot(N,Z1[i]))/c))
            Z2.append(Z1[i] + (g-1)*np.dot(Z1[i],N)*N - g*(a1[i]*v*N)/c)

        px,py,pz = np.array(Z2).T

        w = self.w
        self.update(w,x,y,z,px,py,pz,t)

    def _discretize(self,with_time,queue=None,**kargs):
        """
        See discretize

        Parameters
        ----------
        with_time : bool
            discretize with time
        queue : multiprocessing.Queue
            queue to put results. If none, no multiprocessing is used
        """
        hn=self._ps.hist.hn

        # Initialize
        w       = []
        x,y,z   = [],[],[]
        px,py,pz= [],[],[]
        t       = []

        # Get histo and bins
        if with_time:
            bi,hi=hn([self.x,self.y,self.z,self.px,self.py,self.pz,self.t],normed=False,**kargs)
            bx,by,bz,bpx,bpy,bpz,bt=bi
        else:
            bi,hi=hn([self.x,self.y,self.z,self.px,self.py,self.pz],normed=False,**kargs)
            bx,by,bz,bpx,bpy,bpz=bi

        # Get configurations with non-zero weight
        hb=hi.nonzero()

        # Get data with non-zero weight
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

        if queue is None:
            return w,x,y,z,px,py,pz,t
        else:
            queue.put([w,x,y,z,px,py,pz,t])

    def discretize(self,with_time=True,split=4,MP=True,verbose=True,**kargs):
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
        hist.hn
        """
        hn=self._ps.hist.hn

        if verbose : print('Discretize %s phase space with a bin range splited %i times per axis ...'%(self._ps.specie["name"],split))
        if MP:
            import time
            import multiprocessing as mp
            q = mp.Queue()

        # Initialize
        w       = []
        x,y,z   = [],[],[]
        px,py,pz= [],[],[]
        t       = []

        # Get current brange or reconstruct the default one
        brange = kargs.get('brange',None)
        if brange is None:
            if with_time:
                brange = np.array([
                    [min(self.x),max(self.x)],
                    [min(self.y),max(self.y)],
                    [min(self.z),max(self.z)],
                    [min(self.px),max(self.px)],
                    [min(self.py),max(self.py)],
                    [min(self.pz),max(self.pz)],
                    [min(self.t),max(self.t)]
                ])
            else:
                brange = np.array([
                    [min(self.x),max(self.x)],
                    [min(self.y),max(self.y)],
                    [min(self.z),max(self.z)],
                    [min(self.px),max(self.px)],
                    [min(self.py),max(self.py)],
                    [min(self.pz),max(self.pz)]
                ])
        else:
            del kargs['brange']
        brange = np.array(brange)
        # Calculate width
        width = [(b[1]-b[0])/float(split) for b in brange]
        # Create brange configurations
        brconf=[]
        for ix in range(split):
            # x
            for iy in range(split):
                # y
                for iz in range(split):
                    # z
                    for ipx in range(split):
                        # px
                        for ipy in range(split):
                            # py
                            for ipz in range(split):
                                # pz
                                if with_time:
                                    for it in range(split):
                                        # t
                                        brconf.append([
                                        [brange[0][0] + ix*width[0],brange[0][0] + (ix+1)*width[0]],
                                        [brange[1][0] + iy*width[1],brange[1][0] + (iy+1)*width[1]],
                                        [brange[2][0] + iz*width[2],brange[2][0] + (iz+1)*width[2]],
                                        [brange[3][0] + ipx*width[3],brange[3][0] + (ipx+1)*width[3]],
                                        [brange[4][0] + ipy*width[4],brange[4][0] + (ipy+1)*width[4]],
                                        [brange[5][0] + ipz*width[5],brange[5][0] + (ipz+1)*width[5]],
                                        [brange[6][0] + it*width[6],brange[6][0] + (it+1)*width[6]]
                                        ])
                                else:
                                    brconf.append([
                                    [brange[0][0] + ix*width[0],brange[0][0] + (ix+1)*width[0]],
                                    [brange[1][0] + iy*width[1],brange[1][0] + (iy+1)*width[1]],
                                    [brange[2][0] + iz*width[2],brange[2][0] + (iz+1)*width[2]],
                                    [brange[3][0] + ipx*width[3],brange[3][0] + (ipx+1)*width[3]],
                                    [brange[4][0] + ipy*width[4],brange[4][0] + (ipy+1)*width[4]],
                                    [brange[5][0] + ipz*width[5],brange[5][0] + (ipz+1)*width[5]]
                                    ])
        if not MP:
            i = 0
            for brc in brconf:
                # Print advance approximately 50 times in all the process
                if i%(len(brconf)/50)==0: print("%s/%i configurations done ..."%(i,len(brconf)))
                i+=1
                W,X,Y,Z,Px,Py,Pz,T=self._discretize(with_time=with_time,brange=brc,**kargs)
                # Append current phase space configurations to last ones
                w += list(W)
                x += list(X)        ; y += list(Y)      ; z += list(Z)
                px += list(Px)      ; py += list(Py)    ; pz += list(Pz)
                t += list(T)
        else:
            processes = []
            for brc in brconf:
                proc = mp.Process(target=self._discretize, kwargs = dict(with_time=with_time,queue=q,brange=brc,**kargs))
                processes.append(proc)

            # n = 4
            # for i,proc in enumerate(processes):
            #     if i%n==0 and i!=0:
            #         for j in range(1,n+1):
            #             print("Start process %i"%(i-j))
            #             processes[i-j].start()
            #
            #         for j in range(1,n+1):
            #             processes[i-j].join()
            #             print("Process %i finished"%(i-j))
            n = 4
            for i,proc in enumerate(processes):
                for j in range(1,n+1):
                    print("Start process %i"%((i+n)-j))
                    processes[(i+n)-j].start()

                for j in range(1,n+1):
                    processes[(i+n)-j].join()
                    print("Process %i finished"%((i+n)-j))

            for i in range(len(processes)%n):
                j = len(processes)/n
                k = i + j
                print(k)
                processes[k].start()
            for i in range(len(processes)%n):
                j = len(processes)/n
                k = i + j
                print(k)
                processes[k].join()

            # i=0
            # while i<len(processes):
            #     process[i].start()
            #     i+=1
            #     if i%4==0:
            #         process[]
            # n = 4
            # for i in range(len(processes)/n):
            #     for j in range(1,n):
            #         k = i*j
            #         processes[k].start()
            #         print('k=%i'%k)
            #
            #     for j in range(1,n):
            #         k = i*j
            #         processes[k].join()
            #
            # for i in range(len(processes)%n):
            #     k = len(processes)/n + i
            #     processes[k].start()
            # for i in range(len(processes)%n):
            #     k = len(processes)/n + i
            #     processes[k].join()
            #
            # for proc in processes:
            #     proc.start()
            #
            # i=0
            # nproc = len(processes)
            # while i<nproc:
            #     qsize=int(q.qsize())
            #     if qsize>i:
            #         i+=1
            #         print("{:.0F} % done ...".format(float(i)/nproc * 100))
            #     else:
            #         time.sleep(1.)
            #
            # for proc in processes:
            #     proc.join()

            w,x,y,z,px,py,pz,t=np.array([q.get() for p in processes]).T

        if verbose : print('Done !')

        self.update(w,x,y,z,px,py,pz,t,verbose=verbose)
