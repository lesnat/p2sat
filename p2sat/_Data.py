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

    All the attributes can not be overwriden as they are defined as properties.
    Please call the `update` method to update particle phase-space data.
    """
    def __init__(self,PhaseSpace):
        self._ps= PhaseSpace
        self.raw = _Raw(PhaseSpace)
        self.update(None,None,None,None,None,None,None,None,verbose=False)

    def get_ps(self):
        """
        Return current phase space raw data.
        """
        r = self.raw
        return (r.w,r.x,r.y,r.z,r.px,r.py,r.pz,r.t)

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
            ax = eval("self.raw.%s"%axis)
            # Filter the data if needed
            if select is not None:
                ax = self.select(ax,faxes=list(select.keys()),frange=list(select.values()))
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

        if dct["law"] == "mono":
            # Mono
            K = dct["K"]
            axis = np.array([K] * Nconf)
        elif dct["law"] == "uni":
            # Uniform law
            mini = dct["min"]
            maxi = dct["max"]
            axis = np.random.uniform(mini,maxi,Nconf)
        elif dct["law"] == "exp":
            # Exponential law
            K = dct["K"]
            axis = np.random.exponential(K,Nconf)
        elif dct["law"] == "iso":
            # Isotropic law. Only for angles
            axis = self._generate(dict(law="uni",min=0,max=180),Nconf=Nconf,radians=True)
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

    def generate(self,Nconf,Npart,ekin,theta=None,phi=None,x=None,y=None,z=None,r=None,t=None,verbose=True):
        """
        Generate a particle phase space from given laws.

        TODO : UPDATE DOCUMENTATION

        Parameters
        ----------
        Nconf : int
            total number of configurations
        Npart : float
            total number of particles
        ekin : dict
            parameters to generate kinetic energy
        theta : dict, optional
            parameters to generate theta angle distribution. Default is 'iso'
        phi : dict, optional
            parameters to generate phi angle distribution. Default is 'iso'
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
        - 'uni', for a uniform law in a given range. Keywords 'min' and 'max' MUST be given to specify a minimum and maximum
        - 'exp', for exponential distribution. Charasteristic length/time must be given as a value of keyword `x0`/`y0`/`z0`/`t0`
        - 'gauss', for a gaussian distribution. Center of the distribution must be given with keyword 'mu', and standard deviantion with keyword 'sigma'
        - 'grid', for uniformly placed on a given grid. Keywords 'bins' OR 'min' + 'max' + 'Nbins' MUST be given.

        For dict `r`, available laws are :

        - 'uni', for a uniform law in a given range. Keywords 'min' and 'max' MUST be given to specify a minimum and maximum
        - 'gauss', for a gaussian distribution. Center of the distribution must be given with keyword 'mu', and standard deviantion with keyword 'sigma'

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
        With a `PhaseSpace` object instanciated as `eps`, you can generate
        a mono-energetic source in isotropic direction for 1e12 particles represented
        by 1e6 configurations as follows

        >>> eps.data.generate(Nconf=1e6,Npart=1e12,
        ...                  ekin={"law":"mono","K":20.0},
        ...                  theta={"law":"iso"},
        ...                  phi={"law":"iso"})
        ...
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

        # Ensure that Nconf is of type int (for values such as 1e6)
        Nconf = int(Nconf)

        # Generate weights
        weight  = float(Npart)/Nconf
        g_w     = self._generate(dict(law="mono",K=weight),Nconf)

        # Generate energy
        g_ekin  = self._generate(ekin,Nconf)

        # Generate theta angle
        if theta is not None:
            g_theta = self._generate(theta,Nconf,radians=True)
        else:
            g_theta = self._generate(dict(law="iso"),Nconf,radians=True)

        # Generate phi angle
        if phi is not None:
            g_phi   = self._generate(phi,Nconf,radians=True)
        else:
            g_phi = self._generate(dict(law="iso"),Nconf,radians=True) # FIXME: between 0 and pi ?

        # Reconstruct momentum from energy and angle distributions
        mass    = self._ps.particle["mass"]
        g_p     = np.sqrt(g_ekin**2 + 2*g_ekin*mass)
        g_px    = g_p * np.cos(g_theta)
        g_py    = np.sign(g_phi)*np.sqrt((g_p**2 - g_px**2)/(1. + np.tan(g_phi)**2))
        g_pz    = g_py*np.tan(g_phi)

        # Generate positions and time
        if x is not None:
            g_x     = self._generate(x,Nconf)
        else:
            g_x     = self._generate(dict(law="mono",K=0),Nconf)

        if y is not None:
            g_y     = self._generate(y,Nconf)
        else:
            g_y     = self._generate(dict(law="mono",K=0),Nconf)

        if z is not None:
            g_z     = self._generate(z,Nconf)
        else:
            g_z     = self._generate(dict(law="mono",K=0),Nconf)

        if t is not None:
            g_t     = self._generate(t,Nconf)
        else:
            g_t     = self._generate(dict(law="mono",K=0),Nconf)

        # If transverse position is defined, it replaces g_y,g_z
        if r is not None:
            g_r = self._generate(r,Nconf)
            angle = np.random.uniform(0.,2*np.pi,Nconf)
            g_z = g_r * np.cos(angle)
            g_y = g_r * np.sin(angle)

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

        >>> w,x = eps.data.select('w',['x'],[3])

        or filter by a range (Select all the :math:`\\theta` with :math:`E_{kin} \in [0.511,+\infty]` MeV)

        >>> theta,ekin = eps.data.select('theta',['ekin'],[[0.511,None]])

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
        if verbose: print("Filtering %s phase space with axis %s ..."%(self._ps.particle["name"],faxes))
        data = []
        for ax in self.get_ps():
            # Make a copy of faxes and frange for each axes (because they are modified in select)
            data.append(self.select(ax,list(faxes),list(frange),fpp=fpp))

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

    # def lorentz(self,beta_CM,verbose=True):
    #     """
    #     Lorentz-transformate the particle phase-space with given speed of the center of mass.
    #
    #     TODO
    #
    #     References
    #     ----------
    #     https://en.wikipedia.org/wiki/Lorentz_transformation#Transformation_of_other_quantities
    #     """
    #     if verbose: print("Lorentz-transform %s phase-space with center of mass frame moving at %s c."%(self._ps.particle,beta_CM))
    #     # lowercase : scalar, caption : vector
    #
    #     B = -np.array(beta_CM)
    #     b = np.dot(B,B)
    #     N = B/b
    #     c = 2.99792458e8 * 1e6/1e15 # speed of light in um/fs
    #     v = b * c
    #     g = 1./np.sqrt(1-b**2)
    #
    #     # 4 position
    #     a1 = c*r.t
    #     Z1 = np.array([r.x,r.y,r.z]).T
    #
    #     a2,Z2 =[],[]
    #     for i in range(len(a1)):
    #         a2.append(g*(a1[i] - (v*np.dot(N,Z1[i]))/c))
    #         Z2.append(Z1[i] + (g-1)*np.dot(Z1[i],N)*N - g*(a1[i]*v*N)/c)
    #
    #     t = np.array(a2)/c
    #     x,y,z = np.array(Z2).T
    #
    #     # 4 momentum
    #     a1 = self.ekin/c
    #     Z1 = np.array([r.px,r.py,r.pz]).T
    #
    #     a2,Z2 =[],[]
    #     for i in range(len(a1)):
    #         a2.append(g*(a1[i] - (v*np.dot(N,Z1[i]))/c))
    #         Z2.append(Z1[i] + (g-1)*np.dot(Z1[i],N)*N - g*(a1[i]*v*N)/c)
    #
    #     px,py,pz = np.array(Z2).T
    #
    #     w = self.raw.w
    #     self.update(w,x,y,z,px,py,pz,t)

    def _discretize(self,with_time,**kargs):
        """
        See discretize

        Parameters
        ----------
        with_time : bool
            discretize with time
        """
        hn=self._ps.hist.hn
        r = self.raw

        # Initialize
        w       = []
        x,y,z   = [],[],[]
        px,py,pz= [],[],[]
        t       = []

        # Get histo and bins
        if with_time:
            bi,hi=hn([r.x,r.y,r.z,r.px,r.py,r.pz,r.t],normed=False,**kargs)
            bx,by,bz,bpx,bpy,bpz,bt=bi
        else:
            bi,hi=hn([r.x,r.y,r.z,r.px,r.py,r.pz],normed=False,**kargs)
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

        return w,x,y,z,px,py,pz,t

    def _construct_configuration(self,brange,width,id):
        """
        """
        brc = [[brange[j][0] + id[j]*width[j],brange[j][0] + (id[j]+1)*width[j]] for j in range(len(id))]
        return brc

    def discretize(self,with_time=True,split=4,verbose=True,**kargs):
        """
        Discretize the particles phase space in a 6 or 7 D histogram.

        Parameters
        ----------
        with_time : bool, optional
            discretize with time (7D). Default is True
        split : int
            split phase space into `split` parts per axis, to mitigate MemoryError
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
        r = self.raw
        if verbose : print('Discretize %s phase space with a bin range splited %i times per axis ...'%(self._ps.particle["name"],split))

        # Initialize
        w       = []
        x,y,z   = [],[],[]
        px,py,pz= [],[],[]
        t       = []

        # Get current brange or reconstruct the default one
        brange = kargs.get('brange',None)
        if brange is None:
            brange = np.array([[min(ax),max(ax)] for ax in self.get_ps()])
        else:
            del kargs['brange']
        brange = np.array(brange)
        # Calculate width
        width = [(b[1]-b[0])/float(split) for b in brange]
        # Create brange configurations
        brconf=[]
        for ix in range(split):
            # x
            print("Constructing bin range list : %i / %i ..."%(ix,split))
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
                                        id = [ix,iy,iz,ipx,ipy,ipz,it]
                                        brc=self._construct_configuration(brange,width,id)
                                        W,X,Y,Z,Px,Py,Pz,T=self._discretize(with_time=with_time,brange=brc,**kargs)
                                        # Print advance approximately 50 times in all the process
                                        base = [split**k for k in range(7)]
                                        print("%s/%i configurations done ..."%(np.dot(list(reversed(id)),base),split**7))
                                        # Append current phase space configurations to last ones
                                        w += list(W)
                                        x += list(X)        ; y += list(Y)      ; z += list(Z)
                                        px += list(Px)      ; py += list(Py)    ; pz += list(Pz)
                                        t += list(T)
                                else:
                                    id = [ix,iy,iz,ipx,ipy,ipz,None]
                                    brc=self._construct_configuration(brange,width,id)
                                    W,X,Y,Z,Px,Py,Pz,T=self._discretize(with_time=with_time,brange=brc,**kargs)
                                    # Print advance approximately 50 times in all the process
                                    base = [split**k for k in range(6)]
                                    print("%s/%i configurations done ..."%(np.dot(list(reversed(id)),base),split**6))
                                    # Append current phase space configurations to last ones
                                    w += list(W)
                                    x += list(X)        ; y += list(Y)      ; z += list(Z)
                                    px += list(Px)      ; py += list(Py)    ; pz += list(Pz)
                                    t += [0.]*len(W)

        if verbose : print('Done !')

        self.update(w,x,y,z,px,py,pz,t,verbose=verbose)

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
        if verbose:print("Reshaping phase space ...")
        # Get current phase space configurations
        current_w,*current_ps = self.get_ps()
        current_confs = list(zip(*current_ps))

        # Construct a list of all possible configurations
        reshaped_confs = []
        i=0
        Nconfs = len(current_confs)
        for conf in current_confs:
            if verbose and i%(Nconfs/10)==0:print("Constructing configuration %i/%i ..."%(i,Nconfs))
            i+=1
            if conf not in reshaped_confs:
                reshaped_confs.append(conf)

        # Reconstruct the weights of each new configuration
        reshaped_w = np.zeros(len(reshaped_confs))
        # Loop over all the old configurations, and add their weight to new conf
        for i,conf in enumerate(current_confs):
            if verbose and i%(Nconfs/10)==0: print("Reshaping weight %i/%i ..."%(i,Nconfs))
            # OK to get the first index because each new configuration should be unique
            reshaped_id = reshaped_confs.index(conf)
            reshaped_w[reshaped_id] += current_w[i]

        if verbose:print("Done !")

        reshaped_ps = np.transpose(reshaped_confs)
        self.update(reshaped_w,*reshaped_ps,verbose=verbose)

    def rebin_axis(self,axis,brange,bwidth):
        """
        ...
        """
        axis = self.get_axis(axis)
        # Rebin only if there is different values on axis
        if brange[0] == brange[1]:
            rebined_axis = axis
        else:
            # Initialize with nan
            rebined_axis = np.array([np.nan] * len(axis))
            # Create bins array
            bins = np.arange(brange[0],brange[1],bwidth)
            for b in bins:
                id = self.select("id",[axis],[[b,b+bwidth]])
                rebined_axis[id] = np.array([b] * len(id))
        return rebined_axis

    def rebin_ps(self,bwidth=None,reshape=False,verbose=True):
        """
        # rebin, then reshape
        """
        if verbose:print("Rebinning phase space ...")
        # Get current configurations
        w,*ps = self.get_ps()
        # Rebin current configurations
        rebined_ps = []
        for axis in ps:
            brange = [min(axis),max(axis)]
            bwidth = (brange[1]-brange[0])/10
            rebined_ps.append(self.rebin_axis(axis,brange,bwidth))

        if verbose:print("Done !")

        self.update(w,*rebined_ps,verbose=verbose)

        if reshape:
            self.deduplicate_ps(verbose=verbose)
