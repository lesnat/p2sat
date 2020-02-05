#coding:utf8
import numpy as np

class _MetadataPhaseSpace:
    r"""
    """
    def __init__(self, specie, unit_system):
        self._set_unit_system(unit_system)
        self._set_specie_metadata(specie)
        self._set_quantities_metadata()

    def _set_specie_metadata(self, specie):
        r"""
        """
        if specie in ("gamma","g"):
            self.specie = dict(name = "gamma",  label = r"\gamma",  mass = 0.)
        elif specie in ("positron","e+"):
            self.specie = dict(name = "e+",     label = r"e^+",     mass = 511e3)
        elif specie in ("electron","e-"):
            self.specie = dict(name = "e-",     label = r"e^-",     mass = 511e3)
        elif specie in ("muon+","mu+"):
            self.specie = dict(name = "mu+",    label = r"\mu^+",   mass = 105.6e6)
        elif specie in ("muon-","mu-"):
            self.specie = dict(name = "mu-",    label = r"\mu^-",   mass = 105.6e6)
        elif specie in ("proton","p","H+"):
            self.specie = dict(name = "H+",     label = r"H^+",     mass = 938.3e6)
        elif specie in ("neutron","n"):
            self.specie = dict(name = "n",      label = r"n",       mass = 939.6e6)
        elif specie in ("photon","hv","omega"):
            self.specie = dict(name = "photon", label = r"\omega",  mass = 0.)
        else:
            self.specie = dict(name = "unknown",label = r"unknown", mass = 0.)
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("!!!                  Unknown specie.                     !!!")
            print("!!!        Please define specie's mass (in eV)           !!!")
            print("!!!        in dataset.metadata.specie['mass']            !!!")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

    def _set_quantities_metadata(self):
        r"""
        """
        self.quantity = {}
        self.quantity["id"]         = dict(label = r"id",               dimension = "dimless")
        self.quantity["w"]          = dict(label = r"N",                dimension = "dimless")
        self.quantity["x"]          = dict(label = r"x",                dimension = "length")
        self.quantity["y"]          = dict(label = r"y",                dimension = "length")
        self.quantity["z"]          = dict(label = r"z",                dimension = "length")
        self.quantity["px"]         = dict(label = r"p_x",              dimension = "momentum")
        self.quantity["py"]         = dict(label = r"p_y",              dimension = "momentum")
        self.quantity["pz"]         = dict(label = r"p_z",              dimension = "momentum")
        self.quantity["t"]          = dict(label = r"t",                dimension = "time")
        self.quantity["rx"]         = dict(label = r"r_x",              dimension = "length")
        self.quantity["ry"]         = dict(label = r"r_y",              dimension = "length")
        self.quantity["rz"]         = dict(label = r"r_z",              dimension = "length")
        self.quantity["R"]          = dict(label = r"R",                dimension = "length")
        self.quantity["p"]          = dict(label = r"p",                dimension = "momentum")
        self.quantity["thetax"]     = dict(label = r"\theta_x",         dimension = "angle")
        self.quantity["thetay"]     = dict(label = r"\theta_y",         dimension = "angle")
        self.quantity["thetaz"]     = dict(label = r"\theta_z",         dimension = "angle")
        self.quantity["costhetax"]  = dict(label = r"cos(\theta_x)",    dimension = "dimless")
        self.quantity["costhetay"]  = dict(label = r"cos(\theta_y)",    dimension = "dimless")
        self.quantity["costhetaz"]  = dict(label = r"cos(\theta_z)",    dimension = "dimless")
        self.quantity["phix"]       = dict(label = r"\phi_x",           dimension = "angle")
        self.quantity["phiy"]       = dict(label = r"\phi_y",           dimension = "angle")
        self.quantity["phiz"]       = dict(label = r"\phi_z",           dimension = "angle")
        self.quantity["omegax"]     = dict(label = r"\Omega_x",         dimension = "solid_angle")
        self.quantity["omegay"]     = dict(label = r"\Omega_y",         dimension = "solid_angle")
        self.quantity["omegaz"]     = dict(label = r"\Omega_z",         dimension = "solid_angle")
        self.quantity["etot"]       = dict(label = r"E_{tot}",          dimension = "energy")
        self.quantity["ekin"]       = dict(label = r"E_{kin}",          dimension = "energy")
        self.quantity["gamma"]      = dict(label = r"\gamma",           dimension = "dimless")
        self.quantity["gammax"]     = dict(label = r"\gamma_x",         dimension = "dimless")
        self.quantity["gammay"]     = dict(label = r"\gamma_y",         dimension = "dimless")
        self.quantity["gammaz"]     = dict(label = r"\gamma_z",         dimension = "dimless")
        self.quantity["beta"]       = dict(label = r"\beta",            dimension = "dimless")
        self.quantity["betax"]      = dict(label = r"\beta_x",          dimension = "dimless")
        self.quantity["betay"]      = dict(label = r"\beta_y",          dimension = "dimless")
        self.quantity["betaz"]      = dict(label = r"\beta_z",          dimension = "dimless")
        self.quantity["v"]          = dict(label = r"v",                dimension = "length/time")
        self.quantity["vx"]         = dict(label = r"v_x",              dimension = "length/time")
        self.quantity["vy"]         = dict(label = r"v_y",              dimension = "length/time")
        self.quantity["vz"]         = dict(label = r"v_z",              dimension = "length/time")
        self.quantity["ux"]         = dict(label = r"u_x",              dimension = "dimless")
        self.quantity["uy"]         = dict(label = r"u_y",              dimension = "dimless")
        self.quantity["uz"]         = dict(label = r"u_z",              dimension = "dimless")
        self.quantity["m"]          = dict(label = r"m",                dimension = "energy")
        self.quantity["wekin"]      = dict(label = r"(N \times E_{kin})", dimension = "energy")

    def _set_unit_system(self, unit_system):
        r"""
        """
        self.unit = {}
        self.unit["unit_system"] = unit_system
        if unit_system == "SI-eV":
            self.unit["length"]      = dict(label=r"m",         conv=1.)
            self.unit["momentum"]    = dict(label=r"eV/c",      conv=1.)
            self.unit["energy"]      = dict(label=r"eV",        conv=1.)
            self.unit["time"]        = dict(label=r"s",         conv=1.)
            self.unit["angle"]       = dict(label=r"\pi rad",   conv=np.pi)
            self.unit["solid_angle"] = dict(label=r"\pi sr",    conv=np.pi)
            self.unit["length/time"] = dict(label=r"m/s",       conv=1.)
            self.unit["dimless"]     = dict(label="",           conv=1.)
        elif unit_system == "SI-keV":
            self.unit["length"]      = dict(label=r"m",         conv=1.)
            self.unit["momentum"]    = dict(label=r"keV/c",     conv=1e3)
            self.unit["energy"]      = dict(label=r"keV",       conv=1e3)
            self.unit["time"]        = dict(label=r"s",         conv=1.)
            self.unit["angle"]       = dict(label=r"\pi rad",   conv=np.pi)
            self.unit["solid_angle"] = dict(label=r"\pi sr",    conv=np.pi)
            self.unit["length/time"] = dict(label=r"m/s",       conv=1.)
            self.unit["dimless"]     = dict(label="",           conv=1.)
        elif unit_system == "SI-MeV":
            self.unit["length"]      = dict(label=r"m",         conv=1.)
            self.unit["momentum"]    = dict(label=r"MeV/c",     conv=1e6)
            self.unit["energy"]      = dict(label=r"MeV",       conv=1e6)
            self.unit["time"]        = dict(label=r"s",         conv=1.)
            self.unit["angle"]       = dict(label=r"\pi rad",   conv=np.pi)
            self.unit["solid_angle"] = dict(label=r"\pi sr",    conv=np.pi)
            self.unit["length/time"] = dict(label=r"m/s",       conv=1.)
            self.unit["dimless"]     = dict(label="",           conv=1.)
        elif unit_system == "SI-GeV":
            self.unit["length"]      = dict(label=r"m",         conv=1.)
            self.unit["momentum"]    = dict(label=r"GeV/c",     conv=1e9)
            self.unit["energy"]      = dict(label=r"GeV",       conv=1e9)
            self.unit["time"]        = dict(label=r"s",         conv=1.)
            self.unit["angle"]       = dict(label=r"\pi rad",   conv=np.pi)
            self.unit["solid_angle"] = dict(label=r"\pi sr",    conv=np.pi)
            self.unit["length/time"] = dict(label=r"m/s",       conv=1.)
            self.unit["dimless"]     = dict(label="",           conv=1.)
        elif unit_system == "cgs-eV":
            self.unit["length"]      = dict(label=r"cm",        conv=1e-2)
            self.unit["momentum"]    = dict(label=r"eV/c",      conv=1.)
            self.unit["energy"]      = dict(label=r"eV",        conv=1.)
            self.unit["time"]        = dict(label=r"s",         conv=1.)
            self.unit["angle"]       = dict(label=r"\pi rad",   conv=np.pi)
            self.unit["solid_angle"] = dict(label=r"\pi sr",    conv=np.pi)
            self.unit["length/time"] = dict(label=r"cm/s",      conv=1e-2)
            self.unit["dimless"]     = dict(label="",           conv=1.)
        elif unit_system == "cgs-keV":
            self.unit["length"]      = dict(label=r"cm",        conv=1e-2)
            self.unit["momentum"]    = dict(label=r"keV/c",     conv=1e3)
            self.unit["energy"]      = dict(label=r"keV",       conv=1e3)
            self.unit["time"]        = dict(label=r"s",         conv=1.)
            self.unit["angle"]       = dict(label=r"\pi rad",   conv=np.pi)
            self.unit["solid_angle"] = dict(label=r"\pi sr",    conv=np.pi)
            self.unit["length/time"] = dict(label=r"cm/s",      conv=1e-2)
            self.unit["dimless"]     = dict(label="",           conv=1.)
        elif unit_system == "cgs-MeV":
            self.unit["length"]      = dict(label=r"cm",        conv=1e-2)
            self.unit["momentum"]    = dict(label=r"MeV/c",     conv=1e6)
            self.unit["energy"]      = dict(label=r"MeV",       conv=1e6)
            self.unit["time"]        = dict(label=r"s",         conv=1.)
            self.unit["angle"]       = dict(label=r"\pi rad",   conv=np.pi)
            self.unit["solid_angle"] = dict(label=r"\pi sr",    conv=np.pi)
            self.unit["length/time"] = dict(label=r"cm/s",      conv=1e-2)
            self.unit["dimless"]     = dict(label="",           conv=1.)
        elif unit_system == "cgs-GeV":
            self.unit["length"]      = dict(label=r"cm",        conv=1e-2)
            self.unit["momentum"]    = dict(label=r"GeV/c",     conv=1e9)
            self.unit["energy"]      = dict(label=r"GeV",       conv=1e9)
            self.unit["time"]        = dict(label=r"s",         conv=1.)
            self.unit["angle"]       = dict(label=r"\pi rad",   conv=np.pi)
            self.unit["solid_angle"] = dict(label=r"\pi sr",    conv=np.pi)
            self.unit["length/time"] = dict(label=r"cm/s",      conv=1e-2)
            self.unit["dimless"]     = dict(label="",           conv=1.)
        elif unit_system == "UHI":
            self.unit["length"]      = dict(label=r"\mu m",     conv=1e-6)
            self.unit["momentum"]    = dict(label=r"MeV/c",     conv=1e6)
            self.unit["energy"]      = dict(label=r"MeV",       conv=1e6)
            self.unit["time"]        = dict(label=r"ps",        conv=1e-12)
            self.unit["angle"]       = dict(label=r"\pi rad",   conv=np.pi)
            self.unit["solid_angle"] = dict(label=r"\pi sr",    conv=np.pi)
            self.unit["length/time"] = dict(label=r"um/ps",     conv=1e-6/1e-12)
            self.unit["dimless"]     = dict(label="",           conv=1.)
        else:
            raise NameError("Unknow unit system %s "%unit_system)
