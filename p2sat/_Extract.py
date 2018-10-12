#coding:utf8
import numpy as np

class _Extract(object):
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

    def Smilei_Screen_1d(self,Screen,xnorm,wnorm,tnorm,X=0):
        """
        TODO
        """
        w         = []
        x,y,z     = [],[],[]
        px,py,pz  = [],[],[]
        t         = []

        time= Screen().get()['times']

        Px  = Screen().get()['px'] * 0.511
        Py  = Screen().get()['py'] * 0.511

        wNorm = wnorm/(0.511**2)
        wNorm *= (max(Px)-min(Px))/len(Px)
        wNorm *= (max(Py)-min(Py))/len(Py)

        cdata=np.array([[0.]*len(Px)]*len(Py))

        print("Extracting screen data ...")
        for it,et in enumerate(time):
            ldata = cdata
            cdata = Screen(timesteps=et).getData()[0]
            for ipx,epx in enumerate(cdata):
                for ipy,epy in enumerate(epx):
                    depy = epy-ldata[ipx][ipy]
                    if depy!=0.:
                        w.append(depy*wNorm)
                        x.append(X)
                        y.append(0.)
                        z.append(0.)
                        px.append(Px[ipx])
                        py.append(Py[ipy])
                        pz.append(0.)
                        t.append(et*tnorm)


        pz = [0.0] * len(w)
        x = [X] * len(w)
        z = [0.0] * len(w)
        print("Done !")
        self._ps.data.update(w,x,y,z,px,py,pz,t)

    def Geant4_csv(self,file_name,nthreads=1,verbose=True):
        """
        Extract simulation results from a Geant4 NTuple csv output file

        DEPRECATED

        Parameters
        ----------
        file_name : str
            name of the output file. If it ends with '*_t0.*', the number '0' will be replaced by the number of the current thread
        nthreads : int
            total number of threads to consider
        verbose : bool, optional
            verbosity
        """
        raise DeprecationWarning("Deprecated. Use extract.gp3m2_csv")
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

        w   = data[0::8]
        x   = data[1::8]
        y   = data[2::8]
        z   = data[3::8]
        px  = data[4::8]
        py  = data[5::8]
        pz  = data[6::8]
        t   = data[7::8]
        if verbose:print("Done !")

        self._ps.data.update(w,x,y,z,px,py,pz,t,verbose)

    def gp3m2_csv(self,base_name,verbose=True):
        """
        Extract simulation results from a gp3m2 NTuple csv output file

        Parameters
        ----------
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

        >>> eps.extract.gp3m2_csv("Al_target")
        """
        # Get gp3m2 particle name from p2sat particle name
        part = self._ps.particle["name"]
        if part=="e-":
            part_name = "electron"
        elif part=="e+":
            part_name = "positron"
        elif part=="gamma":
            part_name = "gamma"

        # Construct file base name
        fbase = base_name+"_nt_"+part_name+"_t"
        fext = ".csv"

        # Initialize data list
        data = []
        # Loop over threads
        i = 0
        while True:
            fname = fbase + str(i) + fext
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

        Examples
        --------
        >>> eps = p2sat.PhaseSpace(particle="e-")
        >>> eps.extract.TrILEns_output("../TrILEns/")
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
        TODO
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
