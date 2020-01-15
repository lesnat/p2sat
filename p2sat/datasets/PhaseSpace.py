#coding:utf8
import numpy as np

from ._Dataset import _Dataset
from ._EditPhaseSpace import _EditPhaseSpace
from ._ReadPhaseSpace import _ReadPhaseSpace
from ._SavePhaseSpace import _SavePhaseSpace
from ._LoadPhaseSpace import _LoadPhaseSpace


class PhaseSpace(_Dataset):
    """
    Dataset for specie phase-space analysis.

    Parameters
    ----------
    specie : str
        Name of the specie. Availables are gamma,e-,e+,mu-,mu+,proton,neutron.
    label : str, optional
        ...

    Attributes
    ----------
    load : sub-object
        load phase space from a file
    edit : sub-object
        edit the dataset
    read : sub-object
        read the dataset
    save : sub-object
        save phase space into a file

    Examples
    --------
    ...

    Notes
    -----
    See sub-objects documentation for more informations
    """
    def __init__(self, specie=None, unit_system="SI-eV"):
        self.read     = _ReadPhaseSpace(self, specie, unit_system)
        self.edit     = _EditPhaseSpace(self)
        self.load     = _LoadPhaseSpace(self)
        self.save     = _SavePhaseSpace(self)

    def __add__(self,other):
        """
        Return a new PhaseSpace object, combination of the 2 previous.
        """
        part1 = self.read.metadata.specie
        part2 = other.read.metadata.specie
        if part1["name"] != part2["name"]:
            raise NameError("Can't concatenate phase space for different species")

        ps = PhaseSpace(specie=part1["name"])

        rs = self.read
        ro = other.read

        us = self.read.metadata.unit
        uo = other.read.metadata.unit

        w   = np.array([list(rs.w)                          + list(ro.w)])[0]
        x   = np.array([list(rs.x * us["length"]["conv"])   + list(ro.x * uo["length"]["conv"])])[0]
        y   = np.array([list(rs.y * us["length"]["conv"])   + list(ro.y * uo["length"]["conv"])])[0]
        z   = np.array([list(rs.z * us["length"]["conv"])   + list(ro.z * uo["length"]["conv"])])[0]
        px  = np.array([list(rs.px * us["momentum"]["conv"])+ list(ro.px * uo["momentum"]["conv"])])[0]
        py  = np.array([list(rs.py * us["momentum"]["conv"])+ list(ro.py * uo["momentum"]["conv"])])[0]
        pz  = np.array([list(rs.pz * us["momentum"]["conv"])+ list(ro.pz * uo["momentum"]["conv"])])[0]
        t   = np.array([list(rs.t * us["momentum"]["conv"]) + list(ro.t * uo["momentum"]["conv"])])[0]

        ps.edit.update(w,x,y,z,px,py,pz,t, in_code_units=True, verbose=False)

        return ps

    def __str__(self):
        """
        Returns informations about current `PhaseSpace` object.
        """
        txt  = "\n"
        txt += "fp2sat PhaseSpace instance located at %s\n\n"%hex(id(self))
        txt += "Specie                    : %s\n"%self.read.metadata.specie["name"]
        txt += "Number of configurations  : %i\n"%len(self)
        txt += "Total number of particles : %.4E\n\n"%sum(self.read.w)

        txt += "Statistics                : ( min      ,  max      ,  mean     ,  std      ) unit\n"

        for qty_name in self.read.quantity.keys():
            qty_val = self.read.quantity(qty_name)
            qty_dim = self.read.metadata.quantity[qty_name]["dimension"]
            qty_unit = self.read.metadata.unit[qty_dim]
            txt += "    {qty_name} : ({mini: .3E}, {maxi: .3E}, {mean: .3E}, {std: .3E}) {unit}\n".format(
                        qty_name=qty_name.ljust(21),mini=qty_val.min(),maxi=qty_val.max(),mean=qty_val.mean(),std=qty_val.std(),unit=qty_unit)

        return txt

    def copy(self,verbose=False):
        """
        Return a copy of the current PhaseSpace object.

        Parameters
        ----------
        verbose : bool
          verbosity of the function. If True, a message is displayed when the attributes are loaded in memory
        """
        new = PhaseSpace(specie=self.read.metadata.specie["name"], unit_system=self.read.metadata.unit["unit_system"])
        r=self.read
        new.edit.update(r.w,r.x,r.y,r.z,r.px,r.py,r.pz,r.t, in_code_units=False, verbose=verbose)

        return new
