#coding:utf8
import numpy as np
from ._EditCommon import _EditCommon

class _EditPhaseSpace(_EditCommon):
    r"""
    Edit the dataset.

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
        self._ds= PhaseSpace
        self.update(None,None,None,None,None,None,None,None,in_code_units=True,verbose=False)

    def update(self, w, x, y, z, px, py, pz, t, in_code_units, verbose=True):
        r"""
        Update class attributes with new values.

        Parameters
        ----------
        w,x,y,z,px,py,pz,t : list or numpy.ndarray
          particle phase space. More information can be found in data object documentation
        verbose : bool
          verbosity of the function. If True, a message is displayed when the attributes are loaded in memory
        """
        if verbose: print("Updating dataset ...")
        # Ensure the data are numpy arrays
        W   = np.array(w)
        X   = np.array(x)
        Y   = np.array(y)
        Z   = np.array(z)
        Px  = np.array(px)
        Py  = np.array(py)
        Pz  = np.array(pz)
        t   = np.array(t)
        # Convert in code units if needed
        if not in_code_units:
            X  *= self._ds.metadata.unit["length"]["conv"]
            Y  *= self._ds.metadata.unit["length"]["conv"]
            Z  *= self._ds.metadata.unit["length"]["conv"]
            Px *= self._ds.metadata.unit["momentum"]["conv"]
            Py *= self._ds.metadata.unit["momentum"]["conv"]
            Pz *= self._ds.metadata.unit["momentum"]["conv"]
            t  *= self._ds.metadata.unit["time"]["conv"]
        # Save values
        self._ds.read._w  = W
        self._ds.read._x  = X
        self._ds.read._y  = Y
        self._ds.read._z  = Z
        self._ds.read._px = Px
        self._ds.read._py = Py
        self._ds.read._pz = Pz
        self._ds.read._t  = t
        if verbose: print("Done !")

    def generate(self, Nmp, Np, propagation_axis,
                ekin, phi=None, costheta=None,
                x=None, y=None, z=None, t=None,
                seed=None, verbose=True):
        r"""
        Generate a particle phase space from given functions.

        Parameters
        ----------
        Nmp: int
            Number of macro-particles (number of lines in output file).
        Np: float
            Total number of particles to represent (sum of all the weights).
        propagation_axis: str
            Propagation axis. Must be 'x', 'y' or 'z'.
        ekin: function
            Function to call to generate the kinetic energy of macro-particles.
        costheta: function, optional
            Function to call to generate the cosinus of theta angle of macro-particles. Default is 1 (along propagation_axis).
        phi: function, optional
            Function to call to generate the phi angle (in radians) of macro-particles. Default is uniform between 0 and 2*pi (isotropic).
        x, y, z: function, optional
            Function to call to generate the positions of macro-particles. Default is 0.
        t: function, optional
            Function to call to generate the time of macro-particles. Default is 0.

        Notes
        -----
        All the functions must return values in user units.

        Examples
        --------
        >>> import numpy as np
        >>> eps = ExamplePhaseSpace()
        >>> eps.edit.generate(Nmp=1000, Np=1e12, propagation_axis="x", ekin=lambda:np.random.exponential(scale=2.))
        Generating phase space ...
        Done !
        Updating raw values ...
        Done !

        """
        # Print a starting message
        if verbose: print("Generating phase space ...")

        # Set the random seed
        np.random.seed(seed)

        # Ensure that Nmp is of type int (for values such as 1e6)
        Nmp = int(Nmp)

        # Define default behaviour
        if costheta is None: costheta = lambda: 1.
        if phi is None: phi = lambda: np.random.uniform(0, 2*np.pi)
        if x is None: x = lambda: 0.
        if y is None: y = lambda: 0.
        if z is None: z = lambda: 0.
        if t is None: t = lambda: 0.

        # Generate weights
        g_w         = np.array([float(Np)/Nmp] * Nmp)
        # Generate energy and angles
        g_ekin      = np.array([ekin() for _ in range(Nmp)])
        g_costheta  = np.array([costheta() for _ in range(Nmp)])
        g_phi       = np.array([phi() for _ in range(Nmp)])
        # Generate positions and times
        g_x         = np.array([x() for _ in range(Nmp)])
        g_y         = np.array([y() for _ in range(Nmp)])
        g_z         = np.array([z() for _ in range(Nmp)])
        g_t         = np.array([t() for _ in range(Nmp)])

        # Reconstruct momentum from energy and angle distributions
        mass    = self._ds.metadata.specie["mass"] / self._ds.metadata.unit["energy"]["conv"]
        g_p     = np.sqrt(g_ekin**2 + 2*g_ekin*mass)
        # g_px    = g_p * np.cos(g_theta)
        g_px    = g_p * g_costheta
        g_py_sign = 2 * np.random.randint(0, 2, size = Nmp) - 1
        g_py    = g_py_sign * np.sqrt((g_p**2 - g_px**2)/(1. + np.tan(g_phi)**2))
        g_pz    = g_py*np.tan(g_phi)

        if verbose: print("Done !")

        # Update current object
        if propagation_axis == "x":
            # x y z
            self.update(g_w,g_x,g_y,g_z,g_px,g_py,g_pz,g_t, in_code_units=False, verbose=verbose)
        elif propagation_axis == "y":
            # y z x
            self.update(g_w,g_x,g_y,g_z,g_py,g_pz,g_px,g_t, in_code_units=False, verbose=verbose)
        elif propagation_axis == "z":
            # z x y
            self.update(g_w,g_x,g_y,g_z,g_pz,g_px,g_py,g_t, in_code_units=False, verbose=verbose)

    def filter(self, select, verbose=True):
        r"""
        Filter all the phase space with given condition

        Parameters
        ----------
        select : dict
            filtering dictionary
        verbose : bool, optional
            verbosity

        Examples
        --------
        >>> eps = ExamplePhaseSpace()
        >>> eps.edit.filter(select={'x':[-5.,5.],'r':[0,10],'t':[150,None]})
        Filtering e- phase space with axes ['x', 'r', 't'] ...
        Done !
        Updating raw values ...
        Done !

        >>> eps.read.ekin
        array([1.14810454e-01, 1.78725015e+00, 6.30877382e-01, ..., 3.37300348e+00])
        """
        if verbose: print("Filtering dataset with quantities %s ..."%list(select.keys()))
        dataset_filtered = []
        for qty in ["w","x","y","z","px","py","pz","t"]:
            dataset_filtered.append(self._ds.read.quantity(qty, select))

        if verbose: print("Done !")

        self.update(*dataset_filtered,verbose=verbose)

    def translate(self, Tx=0., Ty=0., Tz=0.,verbose=True):
        r"""
        Translate the particle phase space.

        Parameters
        ----------
        Tx,Ty,Yz : float, optional
            Translate (x,y,z) position of Tx,Ty,Tz. Default is (0,0,0).
        verbose : bool, optional
            Verbosity.
        """
        if verbose: print("Translating phase space by (%.2E,%.2E,%.2E) %s"%(Tx,Ty,Tz,self._ds.metadata.unit["length"]["label"]))

        r = self._ds.read
        dataset_translated = [r.w, x + Tx, y + Ty, z + Tz, r.x, r.py, r.pz, r.t]

        self.update(*dataset_translated, in_code_units=False, verbose=verbose)

    def rotate_x(self, angle, verbose=True):
        r"""
        Rotate the particle phase space along axis x.

        Parameters
        ----------
        angle : float, optional
            rotate (x,y,z) and (px,py,pz) of given angle along x axis.
        verbose : bool, optional
            verbosity

        References
        ----------
        https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
        """
        if verbose:
            print("Rotate phase space of %.2E deg along x axis ..."%angle)

        # Define rotation matrix
        a = np.radians(angle)
        Rx = np.array([
            [1              ,0              ,0              ],
            [0              ,np.cos(a)      ,-np.sin(a)     ],
            [0              ,np.sin(a)      ,np.cos(a)      ]])

        # Loop over particles
        W,X,Y,Z,Px,Py,Pz,t=[],[],[],[],[],[],[],[]
        for w, x, y, z, px, py, pz, t in self._ds.read.particles:
            # Save weight and time
            W.append(w)
            T.append(t)

            # Rotate position
            rX = np.dot([x,y,z], Rx)
            X.append(rX[0])
            Y.append(rX[1])
            Z.append(rX[2])

            # Rotate momentum
            rP = np.dot([px,py,pz], Rx)
            Px.append(rP[0])
            Py.append(rP[1])
            Pz.append(rP[2])

        # Update raw data
        self.update(W,X,Y,Z,Px,Py,Pz,T, in_code_units=False, verbose=verbose)

    def rotate_y(self, angle, verbose=True):
        r"""
        Rotate the particle phase space along axis y.

        Parameters
        ----------
        angle : float, optional
            rotate (x,y,z) and (px,py,pz) of given angle along y axis.
        verbose : bool, optional
            verbosity

        References
        ----------
        https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
        """
        if verbose:
            print("Rotate phase space of %.2E deg along y axis ..."%angle)

        # Define rotation matrix
        a = np.radians(angle)
        Ry = np.array([
            [np.cos(a)      ,0              ,np.sin(a)      ],
            [0              ,1              ,0              ],
            [-np.sin(a)     ,0              ,np.cos(a)      ]])

        # Loop over particles
        W,X,Y,Z,Px,Py,Pz,t=[],[],[],[],[],[],[],[]
        for w, x, y, z, px, py, pz, t in self._ds.read.particles:
            # Save weight and time
            W.append(w)
            T.append(t)

            # Rotate position
            rX = np.dot([x,y,z],Ry)
            X.append(rX[0])
            Y.append(rX[1])
            Z.append(rX[2])

            # Rotate momentum
            rP = np.dot([px,py,pz],Ry)
            Px.append(rP[0])
            Py.append(rP[1])
            Pz.append(rP[2])

        # Update raw data
        self.update(W,X,Y,Z,Px,Py,Pz,T, in_code_units=False, verbose=verbose)

    def rotate_z(self, angle, verbose=True):
        r"""
        Rotate the particle phase space along axis z.

        Parameters
        ----------
        angle : float, optional
            rotate (x,y,z) and (px,py,pz) of given angle along z axis.
        verbose : bool, optional
            verbosity

        References
        ----------
        https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
        """
        if verbose:
            print("Rotate phase space of %.2E deg along z axis ..."%angle)

        # Define rotation matrix
        a = np.radians(angle)
        Rz = np.array([
            [np.cos(a)      ,-np.sin(a)     ,0              ],
            [np.sin(a)      ,np.cos(a)      ,0              ],
            [0              ,0              ,1              ]])

        # Loop over particles
        W,X,Y,Z,Px,Py,Pz,t=[],[],[],[],[],[],[],[]
        for w, x, y, z, px, py, pz, t in self._ds.read.particles:
            # Save weight and time
            W.append(w)
            T.append(t)

            # Rotate position
            rX = np.dot([x,y,z],Rz)
            X.append(rX[0])
            Y.append(rX[1])
            Z.append(rX[2])

            # Rotate momentum
            rP = np.dot([px,py,pz],Rz)
            Px.append(rP[0])
            Py.append(rP[1])
            Pz.append(rP[2])

        # Update raw data
        self.update(W,X,Y,Z,Px,Py,Pz,T, in_code_units=False, verbose=verbose)

    def propagate(self,x=None,t=None,update=True,verbose=True):
        r"""
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

        r = self._ds.read

        W = r.w
        Px = r.px
        Py = r.py
        Pz = r.pz

        if t is not None:
            if verbose: print("Propagate %s phase-space to t = %.4E fs."%(self._ds.particle["name"],t))
            T = np.array([t]*len(W))
            DT = T - r.t
            X = r.x + (r.px/r.p)*r.v*DT
            Y = r.y + (r.py/r.p)*r.v*DT
            Z = r.z + (r.pz/r.p)*r.v*DT

        if x is not None:
            if verbose: print("Propagate %s phase-space to x = %.4E um."%(self._ds.particle["name"],x))
            X = np.array([x]*len(W))
            DT = (X - r.x)/r.v
            T = r.t + DT
            Y = r.y + (r.py/r.p)*r.v*DT
            Z = r.z + (r.pz/r.p)*r.v*DT

        if update:
            self.update(W,X,Y,Z,Px,Py,Pz,T,verbose=verbose)
        else:
            return W,X,Y,Z,Px,Py,Pz,T

    def sort(self, verbose=True):
        r"""
        Sort the particles to have increasing number on x, then on y, then on z, ...
        """
        if verbose:print("Sorting macro-particles ...")
        # Get current dataset
        dataset_unsorted = self._ds.read.dataset

        # Define type of each quantity
        dtype = [('w',np.float64),('x',np.float64),('y',np.float64),('z',np.float64),
        ('px',np.float64),('py',np.float64),('pz',np.float64),('t',np.float64)]

        # Construct list of particles
        particles_unsorted = np.array(list(zip(*dataset_unsorted)), dtype=dtype)

        # Sort particles
        particles_sorted = np.sort(particles_unsorted, order=['x','y','z','px','py','pz','t'])

        # Reconstruct dataset and update
        dataset_sorted = [particles_sorted[label] for label in ['w','x','y','z','px','py','pz','t']]

        if verbose:print("Done !")

        self.update(*dataset_sorted, in_code_units=False, verbose=verbose)

    def _rebin_axis(self,qty,nbins=100,brange=None,queue=None,verbose=True):
        r"""
        Rebin the given qty.

        Parameters
        ----------
        qty : str
            qty to rebin
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
        qty_name = qty
        if verbose: print("Rebinning qty %s ..."%qty_name)
        qty = self._ds.read.quantity(qty_name)

        # Define bin range
        if brange is None: brange = [None, None]
        if brange[0] is None: brange[0] = min(qty)
        if brange[1] is None: brange[1] = max(qty)

        # Rebin only if there is different values on qty
        if brange[0] == brange[1]:
            rebined_axis = qty
        else:
            # Initialize with nan
            rebined_axis = np.array([np.nan] * len(qty))
            # Create bins array
            bwidth = (brange[1]-brange[0])/float(nbins)
            bins = np.linspace(brange[0],brange[1]+bwidth,nbins+1)
            # Loop over all bins
            for i in range(nbins):
                # Filter all the particles that have given qty value in the bin
                id = self.filter_axis("id",select={qty_name:[bins[i],bins[i+1]]},fpp=1e-12)
                # Replace filtered values by the value of bins[i]
                rebined_axis[id] = bins[i]

        if verbose: print("Done !")
        # Put result into given queue or return the rebined qty
        if queue is not None:
            queue.put({qty_name:rebined_axis})
        else:
            return rebined_axis

    def rebin(self, Dx, Dy, Dz, Dpx, Dpy, Dpz, Dt, nbins=100,MP=False,brange=None,deduplicate=False,verbose=True):
        r"""
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
        w = self._ds.data.raw.w
        if type(nbins) is int: nbins = [nbins] * 7
        if brange is None: brange = [[None,None] for _ in range(7)]
        # Rebin current configurations
        ps = []
        if MP:
            import multiprocessing as mp
            q = mp.Manager().Queue()
            processes = []
            for i,ax in enumerate(["x","y","z","px","py","pz","t"]):
                proc = mp.Process(target=self.rebin_axis, args=(ax,), kwargs=dict(nbins=nbins[i],queue=q,brange=brange[i],verbose=verbose))
                processes.append(proc)
                proc.start()

            for proc in processes:
                proc.join()

            ps_dict = {}
            while not q.empty():
                ps_dict.update(q.get())

            self.update(w, ps_dict['x'], ps_dict['y'], ps_dict['z'], ps_dict['px'], ps_dict['py'], ps_dict['pz'], ps_dict['t'], verbose=verbose)

        else:
            for i,ax in enumerate(["x","y","z","px","py","pz","t"]):
                ps.append(self.rebin_axis(ax,nbins=nbins[i],brange=brange[i],verbose=verbose))

            if verbose:print("Done !")
            self.update(w, *ps, verbose=verbose)

        if deduplicate:
            self.deduplicate_ps(verbose=verbose)

    def _merge_hist(self):
        r"""
        """
        Nmps = len(confs_old)
        # Merge
        for i, conf_old in enumerate(confs_old):
            # Continue if already summed
            if verbose and i % (Nmps//100) == 0: print(" %i / %i ..."%(i, Nmps))
            if np.isclose(w_old[i],0.):
                continue
            else:
                # Append current conf to confs_new
                w_new.append(w_old[i])
                confs_new.append(conf_old)
                # Loop over next confs if current conf is not the last one
                j=1
                while i + j < Nmps:
                    if np.allclose(conf_old, confs_old[i+j]):
                        # If next conf is equal to current conf, add its weight to last new conf
                        w_new[-1] += w_old[i+j]
                        # Set old conf weight to 0.
                        w_old[i+j] = 0.
                        # Update j to loop over next conf
                        j += 1
                    else:
                        break

    def merge(self,method='hist',verbose=True,**kargs):
        r"""
        Eliminate twin configurations by summing their weights.

        merge ? compress ? merge_particles ?

        Parameters
        ----------
        algo : {'bf','sort'}, optional
            Algorithm to use. For 'sort' the phase space is assumed to be already sorted.
        verbose : bool, optional
            verbosity

        Notes
        -----
        The 'bf' algorithm stands for brute force.
        The 'sort' algorithm sort the values before merging.

        References
        ----------
        https://en.wikipedia.org/wiki/Packing_problems
        https://en.wikipedia.org/wiki/Set_cover_problem
        """
        if verbose:print("Merging macro-particles ...")
        # Get current macro-particles
        dataset_old = self._ds.read.dataset # Python2 compatibility
        w_old, ps_old = dataset_old[0], dataset_old[1:]
        confs_old = list(zip(*ps_old))

        # Initialize deduplicated confs lists
        w_new = []
        confs_new = []
        if algo == 'bf':
            for id_old, conf_old in enumerate(confs_old):
                if verbose and id_old % (len(confs_old)//100) == 0: print(" %i / %i ..."%(id_old,len(confs_old)))
                try:
                    id_new = confs_new.index(conf_old)
                    w_new[id_new] += w_old[id_old]
                except ValueError:
                    confs_new.append(conf_old)
                    w_new.append(w_old[id_old])
        elif algo == 'sort':
            Nmps = len(confs_old)
            # Merge
            for i, conf_old in enumerate(confs_old):
                # Continue if already summed
                if verbose and i % (Nmps//100) == 0: print(" %i / %i ..."%(i, Nmps))
                if np.isclose(w_old[i],0.):
                    continue
                else:
                    # Append current conf to confs_new
                    w_new.append(w_old[i])
                    confs_new.append(conf_old)
                    # Loop over next confs if current conf is not the last one
                    j=1
                    while i + j < Nmps:
                        if np.allclose(conf_old, confs_old[i+j]):
                            # If next conf is equal to current conf, add its weight to last new conf
                            w_new[-1] += w_old[i+j]
                            # Set old conf weight to 0.
                            w_old[i+j] = 0.
                            # Update j to loop over next conf
                            j += 1
                        else:
                            break

        if verbose:print("Done !")

        ps_new = np.transpose(confs_new)
        self.update(w_new,*ps_new,verbose=verbose)
