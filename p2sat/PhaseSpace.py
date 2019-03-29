#coding:utf8
import numpy as np

from ._Data import _Data
from ._Hist import _Hist
from ._Plot import _Plot
from ._Save import _Save
from ._Load import _Load
from ._Stat import _Stat


class PhaseSpace(object):
    """
    Main class for particle phase-space analysis.

    Parameters
    ----------
    particle : str
        Name of the particle. Availables are gamma,e-,e+,mu-,mu+,proton,neutron.

    Attributes
    ----------
    particle : dict
        contains informations about particle specie, such as name, mass and label in TeX format
    data : sub-object
        contains raw data and methods to manipulate it, such as discretization or transformation.
    hist : sub-object
        make histograms from data
    plot : sub-object
        plot histograms
    load : sub-object
        load phase space from a file
    export : sub-object
        export phase space into a file
    stat : sub-object
        make statistics on particle phase space

    Examples
    --------
    Assuming you already imported p2sat, you can create a PhaseSpace object
    for, let say, electrons, as follows

    >>> eps = PhaseSpace(particle="e-")

    You can then import data from a file, using the `txt` method of sub-object `load`

    >>> eps.load.txt("examples/example.csv")
    Extracting e- phase space from examples/example.csv ...
    Done !
    Updating raw values ...
    Done !

    and look at the imported data

    >>> eps.data.raw.w
    array([1.83014e+08, 1.83014e+08, 1.07523e+07, ..., 7.97784e+06])

    or print general informations about your data set

    >>> print(eps)
    <BLANKLINE>
    p2sat PhaseSpace instance located at 0x...
    <BLANKLINE>
    Specie                    : e-
    Number of configurations  : 4488
    Total number of particles : 2.1463E+11
    <BLANKLINE>
    Statistics                : ( min      ,  max      ,  mean     ,  std      ) unit
        beta                  : ( 5.665E-02,  9.991E-01,  9.494E-01,  1.311E-01) None
        ekin                  : ( 8.218E-04,  1.146E+01,  4.302E+00,  2.922E+00) MeV
        ekin_density          : ( 3.505E+02,  1.386E+09,  7.110E+07,  1.456E+08) MeV
        etot                  : ( 5.118E-01,  1.197E+01,  4.813E+00,  2.922E+00) MeV
        gamma                 : ( 1.002E+00,  2.343E+01,  9.419E+00,  5.718E+00) None
        m                     : ( 5.118E-01,  1.197E+01,  4.813E+00,  2.922E+00) MeV
        omega                 : ( 1.934E-06,  1.256E+01,  6.732E-01,  1.920E+00) sr
        p                     : ( 2.899E-02,  1.196E+01,  4.756E+00,  2.969E+00) MeV/c
        phi                   : (-1.796E+02,  1.800E+02,  4.684E-01,  9.932E+01) deg
        px                    : (-8.822E-01,  1.195E+01,  4.642E+00,  3.050E+00) MeV/c
        py                    : (-2.800E+00,  3.971E+00,  6.467E-02,  5.616E-01) MeV/c
        pz                    : (-3.320E+00,  2.210E+00,  1.490E-02,  5.190E-01) MeV/c
        r                     : ( 1.980E-01,  1.250E+03,  4.085E+01,  7.740E+01) \mu m
        t                     : ( 2.553E+02,  6.315E+03,  1.074E+03,  4.206E+02) fs
        theta                 : ( 4.495E-02,  1.767E+02,  1.720E+01,  2.571E+01) deg
        ux                    : (-9.984E-01,  1.000E+00,  8.929E-01,  3.056E-01) None
        uy                    : (-9.892E-01,  9.681E-01,  1.390E-02,  2.409E-01) None
        uz                    : (-9.906E-01,  9.633E-01,  4.227E-03,  2.262E-01) None
        v                     : ( 1.698E-02,  2.995E-01,  2.846E-01,  3.932E-02) \mu m/fs
        w                     : ( 3.806E+04,  2.989E+09,  4.782E+07,  1.682E+08) None
        x                     : (-6.776E-15,  3.000E+02,  2.197E+02,  7.828E+01) \mu m
        y                     : (-7.697E+02,  7.961E+02,  3.563E+00,  6.181E+01) \mu m
        z                     : (-1.025E+03,  6.236E+02,  6.436E-02,  6.187E+01) \mu m
    <BLANKLINE>

    You can also make histograms, plots or statistics ...

    Notes
    -----
    See sub-objects documentation for more informations
    """
    def __init__(self,particle):
        self.particle = {}
        if particle in ("gamma","g"):
            self.particle["name"] = "gamma"
            self.particle["mass"] = 0
            self.particle["label"]= r"\gamma"
        elif particle in ("positron","e+"):
            self.particle["name"] = "e+"
            self.particle["mass"] = 0.511
            self.particle["label"]= r"e^+"
        elif particle in ("electron","e-"):
            self.particle["name"] = "e-"
            self.particle["mass"] = 0.511
            self.particle["label"]= r"e^-"
        elif particle in ("muon+","mu+"):
            self.particle["name"] = "mu+"
            self.particle["mass"] = 105.6
            self.particle["label"]= r"\mu^+"
        elif particle in ("muon-","mu-"):
            self.particle["name"] = "mu-"
            self.particle["mass"] = 105.6
            self.particle["label"]= r"\mu^-"
        elif particle in ("proton","p"):
            self.particle["name"] = "proton"
            self.particle["mass"] = 938.3
            self.particle["label"]= r"p"
        elif particle in ("neutron","n"):
            self.particle["name"] = "neutron"
            self.particle["mass"] = 939.6
            self.particle["label"]= r"n"
        else:
            raise NameError("Unknown particle specie.")
        self.data     = _Data(self)
        self.hist     = _Hist(self)
        self.plot     = _Plot(self)
        self.load     = _Load(self)
        self.save     = _Save(self)
        self.stat     = _Stat(self)

    def __add__(self,other):
        """
        Return a new PhaseSpace object, combination of the 2 previous.
        """
        spec1 = self.particle["name"]
        spec2 = other.particle["name"]
        if spec1 != spec2:
            raise NameError("Can't combine phase space for different species")

        ps = PhaseSpace(particle=spec1)

        rs = self.data.raw
        ro = other.data.raw

        w   = np.array([list(rs.w) + list(ro.w)])[0]
        x   = np.array([list(rs.x) + list(ro.x)])[0]
        y   = np.array([list(rs.y) + list(ro.y)])[0]
        z   = np.array([list(rs.z) + list(ro.z)])[0]
        px  = np.array([list(rs.px) + list(ro.px)])[0]
        py  = np.array([list(rs.py) + list(ro.py)])[0]
        pz  = np.array([list(rs.pz) + list(ro.pz)])[0]
        t   = np.array([list(rs.t) + list(ro.t)])[0]

        ps.data.update(w,x,y,z,px,py,pz,t,verbose=False)

        return ps

    def __str__(self):
        """
        Returns informations about current `PhaseSpace` object.
        """
        txt  = "\n"
        txt += "p2sat PhaseSpace instance located at %s\n\n"%hex(id(self))
        txt += "Specie                    : %s\n"%self.particle["name"]
        txt += "Number of configurations  : %i\n"%len(self)
        txt += "Total number of particles : %.4E\n\n"%sum(self.data.raw.w)

        txt += "Statistics                : ( min      ,  max      ,  mean     ,  std      ) unit\n"

        for axis,unit in sorted(self.data.raw.units.items()):
            ax = self.data.get_axis(axis)
            txt += "    {ax} : ({mini: .3E}, {maxi: .3E}, {mean: .3E}, {std: .3E}) {unit}\n".format(
                        ax=axis.ljust(21),mini=ax.min(),maxi=ax.max(),mean=ax.mean(),std=ax.std(),unit=unit)

        return txt

    def __len__(self):
        """
        Return the total number of configurations.
        """
        return len(self.data.raw.w)

    def copy(self,verbose=False):
        """
        Return a copy of the current PhaseSpace object.

        Parameters
        ----------
        verbose : bool
          verbosity of the function. If True, a message is displayed when the attributes are loaded in memory
        """
        new = PhaseSpace(particle=self.particle["name"])
        r=self.data.raw
        new.data.update(r.w,r.x,r.y,r.z,r.px,r.py,r.pz,r.t,verbose=verbose)

        return new


class ExamplePhaseSpace(PhaseSpace):
    """
    Example phase space.
    """
    def __init__(self):
        super(ExamplePhaseSpace, self).__init__(particle="electron")
        ekin = {"law":"exp","scale":1.0}
        theta = {"law":"gauss","mu":1.0,"sigma":5.0}
        phi = {"law":"uni","min":0.0,"max":360}
        x = {"law":"uni","min":-10.0,"max":10.0}
        r = {"law":"gauss","mu":0.0,"sigma":10.0}
        t = {"law":"exp","scale":100.0}
        self.data.generate(Nconf=500,Npart=1e10,
                            ekin=ekin,theta=theta,phi=phi,
                            x=x,r=r,t=t,
                            seed=113,verbose=False)
