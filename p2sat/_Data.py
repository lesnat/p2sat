#coding:utf8
import numpy as np
from ._Raw import _Raw

class _Data(object):
    """
    Class containing raw data and methods to manipulate it.

    Attributes
    ----------
    raw : sub-object
        class containing raw data and physical quantities calculations

    Notes
    -----
    Units :

    - lengths are defined in :math:`10^{-6}` meters (um)
    - momentums are defined in :math:`10^{6}` electron-volt/speed of light (MeV/c)
    - times are defined in :math:`10^{-15}` seconds (fs)
    - energies are defined in :math:`10^{6}` electron-volt (MeV)
    - angles are defined in degrees (deg)

    As all the calculations are done with the previously defined units,
    the input data might be firstly converted to those units.

    All the attributes can not be overwritten as they are defined as properties.
    Please call the `update` method to update particle phase-space data.
    """
    def __init__(self,PhaseSpace):
        self._ps= PhaseSpace
        self.raw = _Raw(PhaseSpace)
        self.update(None,None,None,None,None,None,None,None,verbose=False)

    def get_ps(self,select=None):
        """
        Return current phase space raw data.
        """
        # Filter the data if needed
        if select is not None:
            ps = self.filter_ps(select=select,update=False)
        else:
            r = self.raw
            ps = (r.w,r.x,r.y,r.z,r.px,r.py,r.pz,r.t)

        return ps

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
        >>> eps = ExamplePhaseSpace()
        >>> eps.data.get_axis("x")
        array([-4.26043957, -8.225104  ,  0.25424565, ..., -3.19180518])

        >>> eps.data.get_axis("theta",select={'x':[0.,1.],'ekin':[0.511,None]})
        array([4.85380308, 3.79207276, 0.23348689, 0.59771946, 1.02382589,
           1.65421016, 6.84567286, 4.75129691, 3.82330291, 4.15075404,
           7.23677374, 2.69007983])
        """
        if type(axis) is not str:
            # If not a str, return axis
            return axis # FIXME : OK or not ?
        else:
            # Evaluate axis
            ax = eval("self.raw.%s"%axis)
            # Filter the data if needed
            if select is not None:
                ax = self.filter_axis(ax,select=select)
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
        self.raw._w  = np.array(w)
        self.raw._x  = np.array(x)
        self.raw._y  = np.array(y)
        self.raw._z  = np.array(z)
        self.raw._px = np.array(px)
        self.raw._py = np.array(py)
        self.raw._pz = np.array(pz)
        self.raw._t  = np.array(t)
        if verbose: print("Done !")

    def _generate(self,dct,Nconf,radians=False):
        """
        Generate axis configuration from given law.
        """
        # Convert dict values into radians
        if radians:
            convval = []
            for val in dct.values():
                if isinstance(val,(float,int)):
                    convval.append(np.radians(val))
                else:
                    convval.append(val)

            dct = dict(zip(dct.keys(),convval))

        if dct["law"] == "dirac":
            # Dirac distribution
            center = dct["center"]
            axis = np.array([center] * Nconf)
        elif dct["law"] == "uni":
            # Uniform law
            mini = dct["min"]
            maxi = dct["max"]
            axis = np.random.uniform(mini,maxi,Nconf)
        elif dct["law"] == "exp":
            # Exponential law
            scale = dct["scale"]
            axis = np.random.exponential(scale,Nconf)
        elif dct["law"] == "gauss":
            # Gaussian law
            mu = dct["mu"]
            sigma = dct["sigma"]
            axis = np.random.normal(mu,sigma,Nconf)
        elif dct["law"] == "grid":
            # Points placed on a grid
            mini = dct["min"]
            maxi = dct["max"]
            Nbins = dct['Nbins']
            bins = np.linspace(mini,maxi,Nbins)
            index = np.random.randint(0,len(bins),Nconf)
            axis = bins[index]
        elif dct["law"] == "user":
            # User definition
            bins = dct["bins"]
            index = np.random.randint(0,len(bins),Nconf)
            axis = bins[index]

        return axis

    def generate(self,Nconf,Npart,
                ekin,phi=None,theta=None,omega=None,x=None,y=None,z=None,r=None,t=None,
                seed=None,verbose=True):
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
        theta : dict, optional
            parameters to generate theta angle distribution. Default is isotropic
        phi : dict, optional
            parameters to generate phi angle distribution. Default is isotropic
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

        Available laws are :

        - 'dirac', for a dirac distribution. The center of the distribution should be given as the value of keyword 'center'.
        - 'uni', for a uniform law in a given range. Keywords 'min' and 'max' must be given to specify a minimum and maximum.
        - 'exp', for an exponential distribution. The scale of the exponential should be given as a value of keyword 'scale'.
        - 'gauss', for a gaussian (normal) distribution. Center of the distribution must be given with keyword 'mu', and standard deviantion with keyword 'sigma'.
        - 'grid', for uniformly placed on a given grid. Keywords 'min', 'max' and 'Nbins' must be given.
        - 'user', for a user-defined distribution. Keyword 'bins' must be given with a list of all available possibilities.

        For theta and phi, the defaults are isotropic laws.
        They are created using theta={'law':'uni','min':0,'max':180}
        and phi={'law':'uni','min':0,'max':360}.


        Details of the calculations :

        Considering :math:`E_T = E_k + m_0` being the total energy, with
        :math:`E_k` the kinetic energy and :math:`m_0` the rest mass energy.

        We also have :math:`E_T^2 = p^2 + m_0^2` and :math:`p^2=p_x^2 + p_y^2 + p_z^2`
        with :math:`p` in MeV/c.

        Assuming :math:`\cos{\\theta}=\\frac{p_x}{p}` and
        :math:`\\tan{\\phi}=\\frac{p_z}{p_y}` we finaly get

        - :math:`p = E_k^2 - 2 E_k m_0`
        - :math:`p_x = p \cos{\\theta}`
        - :math:`p_y = \\frac{\phi}{\mid \phi \mid} \sqrt{\\frac{p^2 - p_x^2}{1 + \\tan^2{\phi}}}`
        - :math:`p_z = p_y \\tan{\phi}`

        Examples
        --------
        You can generate a mono-energetic source at 20 MeV with a gaussian
        spreading of 5 degrees on theta, representing 1e12 particles
        by 1e6 configurations as follows :

        >>> eps = ExamplePhaseSpace()
        >>> eps.data.generate(Nconf=1e6,Npart=1e12,
        ...                  ekin=dict(law="dirac",center=20.0),
        ...                  theta=dict(law="gauss",mu=0,sigma=5),
        ...                  seed=113)
        ...
        Generate e- phase-space ...
            ekin  : {'law': 'dirac', 'center': 20.0}
            theta : {'law': 'gauss', 'mu': 0, 'sigma': 5}
            phi   : None
        Done !
        Updating raw values ...
        Done !

        >>> eps.data.raw.theta
        array([4.40494361, 1.56579907, 3.67853302, ..., 5.42816942])
        """
        # Print a starting message
        if verbose:
            print("Generate %s phase-space ..."%(self._ps.particle["name"]))
            print("    ekin  : %s"%ekin)
            print("    theta : %s"%theta)
            print("    phi   : %s"%phi)
            if x is not None: print("    x     : %s"%x)
            if y is not None: print("    y     : %s"%y)
            if z is not None: print("    z     : %s"%z)
            if r is not None: print("    r     : %s"%r)
            if t is not None: print("    t     : %s"%t)

        # Set the random seed
        np.random.seed(seed)

        # Ensure that Nconf is of type int (for values such as 1e6)
        Nconf = int(Nconf)

        # Generate weights
        weight  = float(Npart)/Nconf
        g_w     = self._generate(dict(law="dirac",center=weight),Nconf)

        # Generate energy
        g_ekin  = self._generate(ekin,Nconf)

        # Generate phi angle
        if phi is not None:
            g_phi   = self._generate(phi,Nconf,radians=True)
        else:
            g_phi = self._generate(dict(law="uni",min=0,max=360),Nconf,radians=True) # FIXME: between 0 and pi ?

        # Generate theta angle
        if theta is not None:
            g_theta = self._generate(theta,Nconf,radians=True)
        else:
            g_theta = self._generate(dict(law="uni",min=0,max=180),Nconf,radians=True)

        # Generate theta from omega
        if omega is not None:
            g_omega = abs(self._generate(omega,Nconf))
            g_theta = np.arccos(1 - g_omega/(2*np.pi))

        # Reconstruct momentum from energy and angle distributions
        mass    = self._ps.particle["mass"]
        g_p     = np.sqrt(g_ekin**2 + 2*g_ekin*mass)
        g_px    = g_p * np.cos(g_theta)
        g_py_sign = 2 * np.random.randint(0, 2, size = Nconf) - 1
        g_py    = g_py_sign * np.sqrt((g_p**2 - g_px**2)/(1. + np.tan(g_phi)**2))
        # g_py    = np.sign(g_phi)*np.sqrt((g_p**2 - g_px**2)/(1. + np.tan(g_phi)**2))
        # g_py_sign = 2 * np.random.randint(0,2) - 1
        g_pz    = g_py*np.tan(g_phi)

        # Generate positions and time
        if x is not None:
            g_x     = self._generate(x,Nconf)
        else:
            g_x     = self._generate(dict(law="dirac",center=0),Nconf)

        if y is not None:
            g_y     = self._generate(y,Nconf)
        else:
            g_y     = self._generate(dict(law="dirac",center=0),Nconf)

        if z is not None:
            g_z     = self._generate(z,Nconf)
        else:
            g_z     = self._generate(dict(law="dirac",center=0),Nconf)

        if t is not None:
            g_t     = self._generate(t,Nconf)
        else:
            g_t     = self._generate(dict(law="dirac",center=0),Nconf)

        # If transverse position is defined, it replaces g_y,g_z
        if r is not None:
            g_r = self._generate(r,Nconf)
            angle = np.random.uniform(0.,2*np.pi,Nconf)
            g_z = g_r * np.cos(angle)
            g_y = g_r * np.sin(angle)

        if verbose: print("Done !")

        # Update current object
        self.update(g_w,g_x,g_y,g_z,g_px,g_py,g_pz,g_t,verbose=verbose)

    def filter_axis(self,axis,select,fpp=1e-7):
        """
        Filter an axis with a value/range on another axis.

        Parameters
        ----------
        axis : str or numpy.ndarray
            axis to filter
        select : dict
            filtering dictionary
        fpp : float, optional
            relative floating point precision. Default is 1e-7

        Returns
        -------
        axis : numpy.ndarray
            filtered axis

        Examples
        --------
        It is possible to filter an axis by a range, like for example
        select all the :math:`\\theta` with :math:`E_{kin} \in [0.511,+\infty]` MeV

        >>> eps = ExamplePhaseSpace()
        >>> ftheta = eps.data.filter_axis('theta',select={'ekin':[0.511,None]})
        >>> ftheta
        array([ 4.42559421,  4.85380308,  0.57311281, ..., 0.10975927])

        or with several ranges, like in a given space-time volume
        >>> fekin = eps.data.filter_axis('ekin',select={'x':[-5.,5.],'r':[0,10],'t':[150,None]})
        >>> fekin
        array([1.14810454e-01, 1.78725015e+00, 6.30877382e-01, ..., 3.37300348e+00])

        It is also possible to filter with an int or a float.

        See Also
        --------
        get_axis()
        """
        # Get a copy of axis
        axis=self.get_axis(axis,select=None)

        # Construct filtering axes and filtering ranges
        faxes = []
        frange= []
        for key,val in select.items():
            fax = self.get_axis(key,select=None) # Current filtering axis
            faxes.append(fax)
            # Convert int or float values into a list
            if type(val) in (int,float):
                frange.append([val,val])
            else:
                fra = list(val) # Current filtering range is a copy of select values
                # Default ranges are min and max of current filtering axis
                if fra[0] is None: fra[0] = min(fax)
                if fra[1] is None: fra[1] = max(fax)
                frange.append(fra)

        # Before filtering, consider all the configurations
        filtr = np.array([True] * len(axis))
        # Loop over filtering axes
        for i,_ in enumerate(faxes):
            # Then do a boolean AND between last filtr and current one
            # A filtr element is True when the element of filtering axis is contained in filtering range
            filtr *= np.array([x>frange[i][0]*(1-fpp) and x<frange[i][1]*(1+fpp) for x in faxes[i]])

        # Finally return the filtered axis
        return axis[filtr]

    def filter_ps(self,select,fpp=1e-7,update=False,verbose=True):
        """
        Filter all the phase space with given condition

        Parameters
        ----------
        select : dict
            filtering dictionary
        fpp : float, optional
            relative floating point precision. Default is 1e-7
        update : bool, optional
            update or not the current `PhaseSpace` instance. Default is False
        verbose : bool, optional
            verbosity

        Examples
        --------
        >>> eps = ExamplePhaseSpace()
        >>> eps.data.filter_ps(select={'x':[-5.,5.],'r':[0,10],'t':[150,None]}, update=True)
        Filtering e- phase space with axes ['x', 'r', 't'] ...
        Done !
        Updating raw values ...
        Done !

        >>> eps.data.raw.ekin
        array([1.14810454e-01, 1.78725015e+00, 6.30877382e-01, ..., 3.37300348e+00])

        See Also
        --------
        filter_axis()
        """
        if verbose: print("Filtering %s phase space with axes %s ..."%(self._ps.particle["name"],list(select.keys())))
        data = []
        for ax in self.get_ps():
            data.append(self.filter_axis(ax,select,fpp=fpp))

        if verbose: print("Done !")

        if update:
            self.update(*data,verbose=verbose)
        else:
            return data

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
            process rotation before translation. Default is False
        verbose : bool, optional
            verbosity
        """
        # Set defaults
        if T is None: T = (0.,0.,0.)
        if R is None: R = (0.,0.,0.)
        if verbose:
            print("Transformate %s phase space ..."%self._ps.particle["name"])
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

    def propagate(self,x=None,t=None,update=True,verbose=True):
        """
        Propagate the phase space to a given position or time.

        Parameters
        ----------
        x : float, optional
            propagate the phase-space to position x. Default is None (no propagation)
        t : float, optional
            propagate the phase-space to time t. Default is None (no propagation)
        verbose : bool, optional
            verbosity

        Notes
        -----
        x and t can not be defined simultaneously.
        """
        if t is None and x is None:
            raise ValueError("You must define t or x.")
        if t is not None and x is not None:
            raise ValueError("t and x can not be defined simultaneously.")

        r = self.raw

        W = r.w
        Px = r.px
        Py = r.py
        Pz = r.pz

        if t is not None:
            if verbose: print("Propagate %s phase-space to t = %.4E fs."%(self._ps.particle["name"],t))
            T = np.array([t]*len(W))
            DT = T - r.t
            X = r.x + (r.px/r.p)*r.v*DT
            Y = r.y + (r.py/r.p)*r.v*DT
            Z = r.z + (r.pz/r.p)*r.v*DT

        if x is not None:
            if verbose: print("Propagate %s phase-space to x = %.4E um."%(self._ps.particle["name"],x))
            X = np.array([x]*len(W))
            DT = (X - r.x)/r.v
            T = r.t + DT
            Y = r.y + (r.py/r.p)*r.v*DT
            Z = r.z + (r.pz/r.p)*r.v*DT

        if update:
            self.update(W,X,Y,Z,Px,Py,Pz,T,verbose=verbose)
        else:
            return W,X,Y,Z,Px,Py,Pz,T

    def rescale_axis(self,axis,scale):
        """
        Multiply given axis by given scale.

        Parameters
        ----------
        axis : str
            axis to rescale
        scale : float
            scaling coefficient
        """
        if verbose:print("Rescaling axis %s ..."%axis)
        axis = self.get_axis(axis)
        ps = []
        for ax in self.get_ps():
            if ax is axis:
                ps.append(scale*ax)
            else:
                ps.append(ax)

        if verbose: print("Done !")
        self.update(*ps,verbose=verbose)

    def round_axis(self,axis,decimals=8,verbose=True):
        """
        Round the given axis.

        Parameters
        ----------
        axis : str
            axis to round
        decimals : int, optional
            number of decimals
        verbose : bool, optional
            verbosity

        Examples
        --------
        >>> eps = ExamplePhaseSpace()
        >>> eps.data.raw.x
        array([-4.26043957, -8.225104  ,  0.25424565, ..., -3.19180518])

        >>> eps.data.round_axis("x",decimals=1)
        Rounding axis x ...
        Done !

        >>> eps.data.raw.x
        array([-4.3, -8.2,  0.3, ..., -3.2])
        """
        if verbose:print("Rounding axis %s ..."%axis)
        axis = self.get_axis(axis)
        ps = []
        for ax in self.get_ps():
            if ax is axis:
                ps.append(np.around(axis,decimals))
            else:
                ps.append(ax)

        if verbose: print("Done !")
        self.update(*ps,verbose=verbose)

    def deduplicate_ps(self,verbose=True):
        """
        Eliminate duplicated configurations by summing their weights.

        Parameters
        ----------
        verbose : bool, optional
            verbosity
        """
        if verbose:print("Deduplicating phase space configurations ...")
        # Get current phase space configurations
        current_data = self.get_ps() # Python2 compatibility
        current_w, current_ps = current_data[0], current_data[1:]
        current_confs = list(zip(*current_ps))

        # Construct a list of all possible configurations
        deduped_confs = []
        i=0
        Nconfs = len(current_confs)
        for conf in current_confs:
            if verbose and i % (Nconfs//10) == 0:
                print("Constructing configuration number %i/%i ..."%(i,Nconfs))
            i+=1
            if conf not in deduped_confs:
                deduped_confs.append(conf)

        # Reconstruct the weights of each new configuration
        deduped_w = np.zeros(len(deduped_confs))
        # Loop over all the old configurations, and add their weight to new conf
        for i,conf in enumerate(current_confs):
            if verbose and i % (Nconfs//10) == 0:
                print("Summing weight %i/%i ..."%(i,Nconfs))
            # OK to get the first index because each new configuration should be unique
            deduped_id = deduped_confs.index(conf)
            deduped_w[deduped_id] += current_w[i]

        if verbose:print("Done !")

        deduped_ps = np.transpose(deduped_confs)
        self.update(deduped_w,*deduped_ps,verbose=verbose)

    def rebin_axis(self,axis,nbins=100):
        """
        Rebin the given axis.

        TODO : why is there still nan after rebinning ?

        Parameters
        ----------
        axis : str
            axis to rebin
        nbins : int
            number of bins to use

        See Also
        --------
        round_axis

        Examples
        --------
        >>> eps = ExamplePhaseSpace()
        >>> eps.data.raw.x
        array([-4.26043957, -8.225104  ,  0.25424565, ..., -3.19180518])

        >>> eps.data.rebin_axis("x")
        >>> eps.data.raw.x
        [1]
        """
        axis = self.get_axis(axis)

        # Define bin range
        brange = [min(axis),max(axis)]

        # Rebin only if there is different values on axis
        if brange[0] == brange[1]:
            rebined_axis = axis
        else:
            # Initialize with nan
            rebined_axis = np.array([np.nan] * len(axis))
            # rebined_axis = axis.copy()
            # Create bins array
            bwidth = (brange[1]-brange[0])/nbins
            bins = np.linspace(brange[0],brange[1],nbins)
            # Loop over all bins
            for b in bins:
                # print(b, b+bwidth)
                # id = self.filter_axis("id",list([axis]),list([[round(b,7),round(b+bwidth,7)]]),fpp=1e-7)
                id = self.filter_axis("id",list([axis]),list([[b,b+bwidth]]))
                rebined_axis[id] = b

        return rebined_axis

    def rebin_ps(self,nbins=100,deduplicate=False,verbose=True):
        """
        Rebin all the phase space

        Parameters
        ----------
        nbins : int or list of 7 int
            number of bins to use
        deduplicate : bool, optional
            call or not the deduplicate_ps function
        verbose : bool, optional
            verbosity

        See Also
        --------
        rebin_axis
        deduplicate_ps

        Examples
        --------
        >>> eps = ExamplePhaseSpace()
        >>> eps.data.rebin_ps(nbins=10,deduplicate=True)
        >>> eps.data.raw.w
        [...]

        """
        if verbose:print("Rebinning phase space ...")
        # Get current configurations
        data = self.get_ps() # Python2 compatibility
        w, ps = data[0], data[1:]
        # Create nbins array
        if type(nbins) is int: nbins = [nbins] * 7
        # Rebin current configurations
        rebined_ps = []
        for i,axis in enumerate(ps):
            brange = [min(axis),max(axis)]
            bwidth = (brange[1]-brange[0])/10
            rebined_ps.append(self.rebin_axis(axis,nbins=nbins[i]))

        if verbose:print("Done !")

        self.update(w,*rebined_ps,verbose=verbose)

        if deduplicate:
            self.deduplicate_ps(verbose=verbose)
