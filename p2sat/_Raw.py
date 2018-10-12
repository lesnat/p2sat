#coding:utf8
import numpy as np

class _Raw(object):
    """
    Class containing raw data and physical quantities calculations.

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
    others : numpy.ndarray
        see also documentation to look at all the quantities calculated from phase-space data

    Notes
    -----
    To add a new quantity, please add a new function to this file with the name
    of the quantity, only the `self` parameter and with the decorator `@property`.
    Please also add label and unit definition of this new quantity in `__init__`.
    """
    def __init__(self,PhaseSpace):
        self._ps= PhaseSpace
        part_label = self._ps.particle['label']
        if self._ps.particle['name']=="gamma":
            etot_label = 'E_{\gamma}'
            ekin_label = 'E_{\gamma}'
        else:
            etot_label = 'E_{Tot}'
            ekin_label = 'E_{kin}'

        self.labels = {}                ; self.units = {}
        l = self.labels                 ; u = self.units

        l['w'] = 'N_{%s}'%part_label    ; u['w'] = None
        l['x'] = 'x'                    ; u['x'] = '\mu m'
        l['y'] = 'y'                    ; u['y'] = '\mu m'
        l['z'] = 'z'                    ; u['z'] = '\mu m'
        l['px'] = 'p_x'                 ; u['px'] = 'MeV/c'
        l['py'] = 'p_y'                 ; u['py'] = 'MeV/c'
        l['pz'] = 'p_z'                 ; u['pz'] = 'MeV/c'
        l['t'] = 't'                    ; u['t'] = 'fs'

        l['r'] = 'r'                    ; u['r'] = '\mu m'
        l['p'] = 'p'                    ; u['p'] = 'MeV/c'

        l['etot'] = etot_label          ; u['etot'] = 'MeV'
        l['ekin'] = ekin_label          ; u['ekin'] = 'MeV'
        l['gamma'] = '\gamma'           ; u['gamma'] = None
        l['beta'] = '\\beta'            ; u['beta'] = None
        l['m'] = 'm'                    ; u['m'] = 'MeV'
        l['v'] = 'v'                    ; u['v'] = '\mu m/fs'
        l['ux'] = 'u_x'                 ; u['ux'] = None
        l['uy'] = 'u_y'                 ; u['uy'] = None
        l['uz'] = 'u_z'                 ; u['uz'] = None

        l['theta'] = '\\theta'          ; u['theta'] = 'deg'
        l['phi'] = '\\phi'              ; u['phi'] = 'deg'
        # l['omega'] = '\Omega'           ; u['omega'] = 'sr'

        l['ekin_density'] = '(N_{%s} %s)'%(part_label,ekin_label)
        u['ekin_density'] = 'MeV'

    @property
    def w(self):
        """
        Particle statistical weight
        """
        return self._w

    @property
    def x(self):
        """
        Particle position in propagation direction x (in um)
        """
        return self._x

    @property
    def y(self):
        """
        Particle position in transverse direction y (in um)
        """
        return self._y

    @property
    def z(self):
        """
        Particle position in transverse direction z (in um)
        """
        return self._z

    @property
    def px(self):
        """
        Particle momentum in propagation direction x (in MeV/c)
        """
        return self._px

    @property
    def py(self):
        """
        Particle momentum in transverse direction y (in MeV/c)
        """
        return self._py

    @property
    def pz(self):
        """
        Particle momentum in transverse direction z (in MeV/c)
        """
        return self._pz

    @property
    def t(self):
        """
        Particle time (in fs)
        """
        return self._t

    @property
    def r(self):
        """
        Particle distance to the propagation direction x (in um)

        Notes
        -----
        r is calculated as follow :

        .. math::
            r = \sqrt{y^2+z^2}
        """
        return np.sqrt(self.y**2+self.z**2)

    @property
    def p(self):
        """
        Particle absolute momentum (in MeV/c)

        Notes
        -----
        p is calculated as follow :

        .. math::
            p = \sqrt{p_x^2+p_y^2+p_z^2}
        """
        return np.sqrt(self.px**2+self.py**2+self.pz**2)

    @property
    def theta(self):
        """
        Particle polar angle (between px and the plane (py,pz)) (in deg)

        Notes
        -----
        theta is calculated as follow :

        .. math::
            \\theta = \\arccos{p_x/p}

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->py, y->pz, z->px).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system
        """
        return np.degrees(np.arccos(self.px/self.p))

    @property
    def phi(self):
        """
        Particle azimutal angle (between py and pz) (in deg)

        Notes
        -----
        phi is calculated as follow :

        .. math::
            \phi = \\arctan{p_z/p_y}

        with using the arctan2 function to have a result between -180 and 180 deg

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->py, y->pz, z->px).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system

        https://en.wikipedia.org/wiki/Atan2
        """
        return np.degrees(np.arctan2(self.pz,self.py))

    # @property
    # def omega(self):
    #     """
    #     Particle solid angle (in sr) ????
    #
    #     Notes
    #     -----
    #     omega is calculated as follow :
    #
    #     .. math::
    #         ???
    #
    #     References
    #     ----------
    #     https://en.wikipedia.org/wiki/Solid_angle
    #     """
    #     # return np.pi*self.r**2/self.x**2
    #     theta = np.radians(self.theta)
    #     phi = np.pi + np.radians(self.phi)
    #     return np.sin(theta) * theta * phi

    @property
    def etot(self):
        """
        Particle total energy (in MeV)

        Notes
        -----
        etot is calculated as follow :

        .. math::
            E_{tot} = \sqrt{p^2+m_0^2}

        with :math:`m_0` being the rest mass energy

        References
        ----------
        https://en.wikipedia.org/wiki/Energy-momentum_relation
        """
        mass = self._ps.particle["mass"]
        return np.sqrt(self.p**2 + mass**2)

    @property
    def ekin(self):
        """
        Particle kinetic energy (in MeV)

        Notes
        -----
        ekin is calculated as follow :

        .. math::
            E_{kin} = E_{tot} - m_0

        with :math:`m_0` being the rest mass energy

        References
        ----------
        https://en.wikipedia.org/wiki/Energy-momentum_relation
        """
        mass = self._ps.particle["mass"]
        return self.etot-mass

    @property
    def gamma(self):
        """
        Particle Lorentz factor

        Notes
        -----
        gamma is calculated as follow :

        .. math::
            \gamma = E_{tot}/m_0

        with :math:`m_0` being the rest mass energy

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        mass = self._ps.particle["mass"]
        return self.etot/mass

    @property
    def beta(self):
        """
        Particle normalized velocity

        Notes
        -----
        beta is calculated as follow :

        .. math::
            \\beta = \sqrt{1 - \\frac{1}{\gamma^2}}

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        return np.sqrt(1.-1./self.gamma**2)

    @property
    def v(self):
        """
        Particle absolute velocity (in um/fs)

        Notes
        -----
        v is calculated as follow :

        .. math::
            v = \\beta \\times c

        with :math:`c` being the speed of light

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        c = 2.99792458e8 * 1e6/1e15 # speed of light in um/fs
        return self.beta * c

    @property
    def ux(self):
        """
        Particle direction on x axis

        Notes
        -----
        ux is calculated as follow :

        .. math::
            u_x = \\frac{p_x}{p}
        """
        return self.px/self.p

    @property
    def uy(self):
        """
        Particle direction on y axis

        Notes
        -----
        uy is calculated as follow :

        .. math::
            u_y = \\frac{p_y}{p}
        """
        return self.py/self.p

    @property
    def uz(self):
        """
        Particle direction on z axis

        Notes
        -----
        uz is calculated as follow :

        .. math::
            u_z = \\frac{p_z}{p}
        """
        return self.pz/self.p

    @property
    def m(self):
        """
        Particle relativistic mass (in MeV)

        Notes
        -----
        m is calculated as follow :

        .. math::
            m = \gamma m_0

        with :math:`m_0` being the rest mass energy

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        mass = self._ps.particle["mass"]
        return self.gamma * mass

    @property
    def ekin_density(self):
        """
        Particle energy density (in MeV)

        Notes
        -----
        ekin_density is calculated as follow :

        .. math::
            E_{kin} \\times w
        """
        return self.ekin * self.w

    # @property
    # def brilliance(self):
    #     return photons/s/mm2/mrad2/0.1%BW.
