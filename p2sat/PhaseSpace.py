#coding:utf8
import numpy as np

from ._Data import _Data
from ._Hist import _Hist
from ._Plot import _Plot
from ._Export import _Export
from ._Extract import _Extract
from ._Stat import _Stat


class PhaseSpace(object):
    """
    Main class for particle phase-space analysis.

    Parameters
    ----------
    specie : str
        Name of the particle specie. Availables are gamma,e-,e+,mu-,mu+.

    Attributes
    ----------
    specie : dict
        contains informations about particle specie, such as name, mass and label in TeX format
    data : sub-object
        contains raw data and methods to manipulate it, such as discretization or transformation.
    hist : sub-object
        make histograms from data
    plot : sub-object
        plot histograms
    extract : sub-object
        load phase space from a file
    export : sub-object
        export phase space into a file
    stat : sub-object
        make statistics on particle phase space

    Examples
    --------
    Assuming you already imported p2sat, you can create a PhaseSpace object
    for, let say, electrons, as follows

    >>> eps = p2sat.PhaseSpace(specie="e-")

    You can then import data from a file, using the `txt` method of sub-object `extract`

    >>> eps.extract.txt("example.csv")

    and look at the imported data

    >>> print(eps.data.w)

    or print general informations about your data set

    >>> print(eps)

    You can also make histograms, plots or statistics ...

    Notes
    -----
    See sub-objects documentation for more informations
    """
    def __init__(self,specie):
        self.specie = {}
        if specie in ("gamma","g"):
            self.specie["name"] = "gamma"
            self.specie["mass"] = 0
            self.specie["label"]= "\gamma"
        elif specie in ("positron","e+"):
            self.specie["name"] = "e+"
            self.specie["mass"] = 0.511
            self.specie["label"]= "e^+"
        elif specie in ("electron","e-"):
            self.specie["name"] = "e-"
            self.specie["mass"] = 0.511
            self.specie["label"]= "e^-"
        elif specie in ("muon+","mu+"):
            self.specie["name"] = "mu+"
            self.specie["mass"] = 105.6
            self.specie["label"]= "\mu^+"
        elif specie in ("muon-","mu-"):
            self.specie["name"] = "mu-"
            self.specie["mass"] = 105.6
            self.specie["label"]= "\mu^-"
        else:
            raise NameError("Unknown particle specie.")
        self.data     = _Data(self)
        self.hist     = _Hist(self)
        self.plot     = _Plot(self)
        self.extract  = _Extract(self)
        self.export   = _Export(self)
        self.stat     = _Stat(self)

    def __add__(self,other):
        """
        Return a new PhaseSpace object, combination of the 2 previous.
        """
        spec1 = self.specie["name"]
        spec2 = other.specie["name"]
        if spec1 != spec2:
            raise NameError("Can't combine phase space for different species")

        ps = PhaseSpace(specie=spec1)

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
        txt += "Specie                    : %s\n"%self.specie["name"]
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
        new = PhaseSpace(specie=self.specie["name"])
        r=self.data.raw
        new.data.update(r.w,r.x,r.y,r.z,r.px,r.py,r.pz,r.t,verbose=verbose)

        return new
