#coding:utf8
import numpy as np

class _Load(object):
    """
    Import data from a file.

    Notes
    -----
    If you want to add a method to import data from another code, you must proceed as follow :

    - Add a method to this object, with the name of your code. It must contains the keyword `self` as a first argument (because of object-oriented paradigm), and all the other parameters you need
    - Get the data from your file and put it in lists or numpy arrays, one line describing one particle
    - Call the `update` method of `data` sub-object (access via `self._ps.data.update`)
    - Please write a documentation and share !

    You can copy-paste the `txt` method to have a basic example of file import.
    """
    def __init__(self,PhaseSpace):
        self._ps=PhaseSpace

    def txt(self,file_name,sep=",",verbose=True):
        """
        Load particle phase space from a text file.

        Parameters
        ----------
        file_name : str
            name of the input file
        sep : str
            character used to separate values. Default is ','
        verbose : bool, optional
            verbosity of the function. If True, a message is displayed when the data is imported

        See Also
        --------
        export.txt
        """
        if verbose: print("Extracting %s phase space from %s ..."%(self._ps.particle["name"],file_name))

        # Initialize data lists
        w         = []
        x,y,z     = [],[],[]
        px,py,pz  = [],[],[]
        t         = []

        # Open file
        with open(file_name,'r') as f:
            # Loop over lines
            for line in f.readlines():
                # If current line is not a comment, save data
                if line[0]!="#":
                    data=line.split(sep)
                    w.append(float(data[0]))
                    x.append(float(data[1]))  ; y.append(float(data[2]))  ; z.append(float(data[3]))
                    px.append(float(data[4])) ; py.append(float(data[5])) ; pz.append(float(data[6]))
                    t.append(float(data[7]))

        if verbose: print('Done !')

        # Save data in PhaseSpace object
        self._ps.data.update(w,x,y,z,px,py,pz,t,verbose=verbose)

    def Smilei_Screen_1d(self,path,nb,r,x=0,verbose=True):
        """
        Extract phase space from Smilei 1D Screen diagnostic.

        Parameters
        ----------
        path : str
            path to the simulation folder
        nb : int
            Screen number
        r : float
            typical radius to consider in transverse direction (in um)
        x : float, optional
            diagnostic position
        verbose : bool, optional

        Notes
        -----
        On a 1D Smilei simulation, a typical DiagScreen must be declared as follows
        ::
            DiagScreen(
                shape               = 'plane',
                point               = [xtarget[1] - 5*um],
                vector              = [1.],
                direction           = 'forward',

                deposited_quantity  = 'weight',
                species             = ['e'],
                axes                = [
                     ['px' , pmin     , pmax     , 301],
                     ['py' , -pmax/5  , pmax/5   , 301]
                ],

                every               = every
            )
        """
        if verbose: print("Extracting screen data from %s ..."%path)

        # Import Smilei simulation
        import happi
        S = happi.Open(path,verbose=False)
        nl = S.namelist

        # Define physical constants
        m_e = 9.11e-31
        epsilon_0 = 8.85e-12
        e = 1.6e-19
        c = 2.99792458e8
        epsilon_0 = 8.854187817e-12
        # Smilei's unit in SI
        Wr = nl.Main.reference_angular_frequency_SI
        Tr = 1/Wr
        Lr = c/Wr
        Pr = 0.511 # MeV/c
        # Calculate normalizations
        nc = m_e * epsilon_0 * (Wr/e)**2
        Lx = nl.Main.grid_length[0] * Lr # Use a try/except ?
        vol = Lx * np.pi * (r * 1e-6)**2
        wnorm = nc * vol # Weight normalization : Before -> in Nc/Pr/Pr, After -> in Number/Pr/Pr
        tnorm = 1e-15/Tr
        xnorm = 1e-6/Lr
        # Save diag position
        xdiag = x

        # Initialize phase space lists
        w         = []
        x,y,z     = [],[],[]
        px,py,pz  = [],[],[]
        t         = []

        # Retrieve Screen data
        times = S.Screen(nb).getTimes()
        timesteps= S.Screen(nb).getTimesteps()

        Px  = S.Screen(nb).getAxis("px") * Pr
        Py  = S.Screen(nb).getAxis("py") * Pr

        # Compensate happi correction on weights
        wnorm /= Pr**2 # Weights are now in Nb/(MeV/c)/(MeV/c) (independant of bin size)
        wnorm *= (max(Px)-min(Px))/len(Px) # Multiply by bin size : weights are now in Nb/(MeV/c)/bin
        wnorm *= (max(Py)-min(Py))/len(Py) # Weight are now in Nb/bin/bin (dependant of bin size, it counts number of particles for given conf)

        # Current data is initialized as an empty matrix
        cdata=np.array([[0.]*len(Px)]*len(Py))

        # Loop over times
        for it,et in enumerate(timesteps):
            ldata = cdata
            # Retrieve data for given time
            cdata = S.Screen(nb,timesteps=et).getData()[0]
            # Loop over px then py
            if verbose and it % (len(times)//10) == 0: print("Retrieving timestep %i/%i ..."%(et,timesteps[-1]))
            for ipx,epx in enumerate(cdata):
                for ipy,epy in enumerate(epx):
                    # Get weight difference for given configuration
                    depy = epy-ldata[ipx][ipy]
                    # If non-zero, save config
                    if depy!=0.:
                        w.append(depy * wnorm)
                        px.append(Px[ipx])
                        py.append(Py[ipy])
                        t.append(times[it] * tnorm)

        # Reconstruct missing data
        pz = [0.0] * len(w)
        x = [xdiag] * len(w)
        y = [0.0] * len(w)
        z = [0.0] * len(w)

        # Update current phase space
        if verbose: print("Done !")
        self._ps.data.update(w,x,y,z,px,py,pz,t)

    def Smilei_TrackParticles(self,path,species,verbose=True):
        """
        Extract phase space from a TrackParticles Smilei diagnostic.

        Parameters
        ----------
        path : str
            path to the simulation folder
        species : str
            name of the specie in the Smilei namelist
        verbose : bool, optional
            verbosity
        """
        if verbose: print("Extracting %s phase space from %s TrackParticles ..."%(self._ps.particle["name"],species))
        # Open simulation
        import happi
        S = happi.Open(path,verbose=False)
        timesteps = S.TrackParticles(species=species,sort=False).getTimesteps()
        dt = S.namelist.Main.timestep

        # Initialize ps list
        w         = []
        x,y,z     = [],[],[]
        px,py,pz  = [],[],[]
        t         = []

        # Get data for each timestep
        for ts in timesteps:
            if verbose:print("Timestep %i/%i ..."%(ts,timesteps[-1]))
            data = S.TrackParticles(species=species,timesteps=ts,sort=False).get()[ts]
            id = data["Id"]
            if len(id[id>0]) == 0: continue
            w += list(data["w"][id>0])
            x += list(data["x"][id>0])
            try:
                y += list(data["y"][id>0])
            except:
                y += [0.]
            try:
                z += list(data["z"][id>0])
            except:
                z += [0.]
            px += list(data["px"][id>0])
            py += list(data["py"][id>0])
            pz += list(data["pz"][id>0])
            t += [ts/dt]*len(id[id>0])

        if verbose: print("Done !")

        self._ps.data.update(w,x,y,z,px,py,pz,t,verbose=verbose)

    def gp3m2_csv(self,path,base_name,verbose=True):
        """
        Extract simulation results from a gp3m2 NTuple csv output file

        Parameters
        ----------
        path : str
            path to the simulation folder
        base_name : str
            base file name
        verbose : bool, optional
            verbosity

        Examples
        --------
        For the gp3m2 output file name `Al_target_nt_electron_t0.csv`, the base_name
        is `Al_target`.
        Assuming a `p2sat.PhaseSpace` object is instanciated for particle `e-` as eps,
        you can import simulation results for all the threads as follows

        >>> eps = ExamplePhaseSpace()
        >>> # eps.extract.gp3m2_csv("Al_target")
        """
        # Get gp3m2 particle name from p2sat particle name
        part = self._ps.particle["name"]
        if part=="e-":
            part_name = "electron"
        elif part=="e+":
            part_name = "positron"
        elif part=="gamma":
            part_name = "gamma"
        elif part=="photon":
            part_name = "OpPhoton"

        # Construct file base name
        fbase = base_name+"_nt_"+part_name+"_t"
        fext = ".csv"

        # Initialize data list
        data = []
        # Loop over threads
        i = 0
        while True:
            fname = path + fbase + str(i) + fext
            i    += 1
            try:
                # Open file for thread i-1
                with open(fname,'r') as f:
                    if verbose:print("Extracting %s ..."%fname)
                    # Loop over lines
                    for line in f.readlines():
                        # Save data if current line is not a comment
                        if line[0]!='#':
                            for e in line.split(','):
                                data.append(float(e))
            # If no more thread, break the loop
            except IOError:
                break

        # Get phase space from data list
        w   = data[0::8]
        x   = data[1::8]
        y   = data[2::8]
        z   = data[3::8]
        px  = data[4::8]
        py  = data[5::8]
        pz  = data[6::8]
        t   = data[7::8]
        if verbose:print("Done !")

        # Save phase space data in PhaseSpace object
        self._ps.data.update(w,x,y,z,px,py,pz,t,verbose)

    def TrILEns_output(self,path,verbose=True):
        """
        Extract simulation results from a TrILEns output.txt file

        Parameters
        ----------
        path : str
            simulation path
        verbose : bool, optional
            verbosity
        """
        particle = self._ps.particle["name"]
        if verbose:print("Extracting {} phase space from {}output.txt ...".format(particle,path))

        # Get TrILEns particle label from p2sat particle name
        if particle == "e-":
            label = "electrons"
        elif particle == "e+":
            label = "positrons"
        elif particle == "gamma":
            label = "photons"

        # Initialize phase space lists
        w         = []
        x,y,z     = [],[],[]
        px,py,pz  = [],[],[]
        t         = []

        # Boolean to extract only the data of correct particle
        is_correct_species=False

        # Open output file
        with open(path+'output.txt','r') as f:
            # 34 first lines are informations about the simulation
            for _ in range(3):
                f.readline()
            line = f.readline()
            if line.split()[1]=="T":
                chi_to_t = True
            else:
                chi_to_t = False
            for _ in range(30):
                f.readline()
            # Loop over data
            for line in f.readlines():
                try:
                    # Photons do not have Chi value
                    if label == "photons":
                        W,X,Y,Z,Px,Py,Pz,Gamma=line.split()
                        chi_to_t = False
                    else:
                        W,X,Y,Z,Px,Py,Pz,Gamma,Chi=line.split()
                    # If correct particle, save data
                    if is_correct_species:
                        w.append(float(W))
                        x.append(float(X))     ; y.append(float(Y))   ; z.append(float(Z))
                        px.append(float(Px)*0.511)   ; py.append(float(Py)*0.511) ; pz.append(float(Pz)*0.511)
                        if chi_to_t:
                            t.append(float(Chi))
                        else:
                            t.append(0.)

                # If current line is a string (not possible to read data), test if particle label in current line
                except ValueError:
                    if label in line.split():
                        is_correct_species = True
                    else:
                        is_correct_species = False

        if verbose:print("Done !")

        # Save data in PhaseSpace object
        self._ps.data.update(w,x,y,z,px,py,pz,t,verbose=verbose)

    def TrILEns_prop_ph(self,path,verbose=True):
        """
        Extract simulation results from a TrILEns prop_ph file

        Parameters
        ----------
        path : str
            simulation path
        verbose : bool, optional
            verbosity
        """
        if self._ps.particle["name"]!="gamma":
            raise NameError("prop_ph.t contains informations about gamma photons ! Current particle name is %s"%self._ps.particle["name"])
        if verbose: print("Extracting %s phase space from %s ..."%(self._ps.particle["name"],path+"prop_ph.t"))

        # Initialize data lists
        w         = []
        x,y,z     = [],[],[]
        px,py,pz  = [],[],[]
        t         = []

        with open(path+"prop_ph.t",'r') as f:
            # First line gives information about time
            line = f.readline()
            if line == "8 1.\n":
                with_time = False
            elif line == "9 1.\n":
                with_time = True
            else:
                raise NameError("Unknown time identifier at line 1 : %s"%line)
            # second line is a legend
            _ = f.readline()
            # Loop over data lines
            for line in f.readlines():
                # If current line is not a comment, save data
                data=line.split()
                w.append(float(data[0]))
                x.append(float(data[1]))  ; y.append(float(data[2]))  ; z.append(float(data[3]))
                px.append(float(data[4])) ; py.append(float(data[5])) ; pz.append(float(data[6]))
                if with_time:
                    t.append(float(data[8]))
                else:
                    t.append(0.)

        if verbose: print('Done !')

        self._ps.data.update(w,x,y,z,px,py,pz,t,verbose=verbose)