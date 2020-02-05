#coding:utf8
import numpy as np

class _ReadPhaseSpace(object):
    r"""
    Read the dataset.

    Attributes
    ----------
    metadata : dict
        ...
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
    def __init__(self, PhaseSpace):
        self._ds      = PhaseSpace

    def quantity(self, qty, select=None):
        r"""
        Return the values of given quantity.

        Parameters
        ----------
        qty : str
            Quantity name.
        select : dict
            Filtering dictionnary.

        Examples
        --------
        >>> eps = ExamplePhaseSpace()
        >>> eps.read.values("x")
        array([-4.26043957, -8.225104  ,  0.25424565, ..., -3.19180518])

        >>> eps.read.values("thetax",select={'x':[0.,1.],'ekin':[0.511,None]})
        array([4.85380308, 3.79207276, 0.23348689, 0.59771946, 1.02382589,
           1.65421016, 6.84567286, 4.75129691, 3.82330291, 4.15075404,
           7.23677374, 2.69007983])
        """
        # Evaluate quantity
        qty = eval("self.%s"%qty)

        # Filter the data if needed
        if select is not None:
            # Construct filtering axes and filtering ranges
            fqties  = []
            franges = []
            # For each key of select dictionnary, get the appropriate filtering quantity
            # For each value of select dictionnary, get the appropriate filtering range
            for key, val in select.items():
                # Get current filtering quantity
                fqty = self.quantity(key, select=None)

                # Construct filtering range of current quantity
                if type(val) in (int,float):
                    # Filtering range with int or float are converted into a list with the two same element
                    frange = [val,val]
                else:
                    # Current filtering range is a copy of select values
                    frange = list(val)
                    # Default ranges : min and max of current filtering quantity
                    if frange[0] is None: frange[0] = min(fqty)
                    if frange[1] is None: frange[1] = max(fqty)

                # Save current filtering quantity and filtering range
                fqties.append(fqty)
                franges.append(frange)

            # Before filtering, all the quantity configurations are considered
            filtr = np.array([True] * len(qty))
            # Loop over filtering axes
            for i,_ in enumerate(fqties):
                # Then do a boolean AND between last filtr and current one
                # A filtr element is True when the element of filtering quantity is contained in filtering range
                filtr *= np.greater_equal(fqties[i], franges[i][0]) # OK with floating point precision ?
                filtr *= np.less_equal(fqties[i], franges[i][1])

            # Filter the quantity
            qty = qty[filtr]

        # return result
        return qty

    @property
    def dataset(self):
        r"""
        ???
        """
        return (self.w, self.x, self.y, self.z, self.px, self.py, self.pz, self.t)

    @property
    def particles(self):
        r"""
        ???
        """
        return list(zip(self.w, self.x, self.y, self.z, self.px, self.py, self.pz, self.t))

    @property
    def id(self):
        r"""
        Particle index.
        """
        return np.array(range(len(self.w)))

    @property
    def w(self):
        r"""
        Particle statistical weight.
        """
        return self._w

    @property
    def x(self):
        r"""
        Particle position in axis x.
        """
        x = self._x / self._ds.metadata.unit["length"]["conv"] # Convert to user unit
        return x

    @property
    def y(self):
        r"""
        Particle position in axis y.
        """
        y = self._y / self._ds.metadata.unit["length"]["conv"] # Convert to user unit
        return y

    @property
    def z(self):
        r"""
        Particle position in axis z.
        """
        z = self._z / self._ds.metadata.unit["length"]["conv"] # Convert to user unit
        return z

    @property
    def px(self):
        r"""
        Particle momentum in direction x.
        """
        px = self._px / self._ds.metadata.unit["momentum"]["conv"] # Convert to user unit
        return px

    @property
    def py(self):
        r"""
        Particle momentum in direction y.
        """
        py = self._py / self._ds.metadata.unit["momentum"]["conv"] # Convert to user unit
        return py

    @property
    def pz(self):
        r"""
        Particle momentum in direction z.
        """
        pz = self._pz / self._ds.metadata.unit["momentum"]["conv"] # Convert to user unit
        return pz

    @property
    def t(self):
        r"""
        Particle time (in fs).
        """
        t = self._t / self._ds.metadata.unit["time"]["conv"] # Convert to user unit
        return t

    @property
    def rx(self):
        r"""
        Particle distance to axis x.

        Notes
        -----
        rx is calculated as follow :

        .. math::
            r_x = \sqrt{y^2+z^2}
        """
        y = self.y * self._ds.metadata.unit["length"]["conv"] # Convert to code unit
        z = self.z * self._ds.metadata.unit["length"]["conv"] # Convert to code unit
        rx = np.sqrt(y**2 + z**2) / self._ds.metadata.unit["length"]["conv"] # Convert to user unit
        return rx

    @property
    def ry(self):
        r"""
        Particle distance to axis y.

        Notes
        -----
        ry is calculated as follow :

        .. math::
            r_y = \sqrt{x^2+z^2}
        """
        x = self.x * self._ds.metadata.unit["length"]["conv"] # Convert to code unit
        z = self.z * self._ds.metadata.unit["length"]["conv"] # Convert to code unit
        ry = np.sqrt(x**2 + z**2) / self._ds.metadata.unit["length"]["conv"] # Convert to user unit
        return ry

    @property
    def rz(self):
        r"""
        Particle distance to axis z.

        Notes
        -----
        rz is calculated as follow :

        .. math::
            r_z = \sqrt{x^2+y^2}
        """
        x = self.x * self._ds.metadata.unit["length"]["conv"] # Convert to code unit
        y = self.y * self._ds.metadata.unit["length"]["conv"] # Convert to code unit
        rz = np.sqrt(x**2 + y**2) / self._ds.metadata.unit["length"]["conv"] # Convert to user unit
        return rz

    @property
    def R(self):
        r"""
        Particle distance to the origin.

        Notes
        -----
        R is calculated as follow :

        .. math::
            R = \sqrt{x^2+y^2+z^2}
        """
        x = self.x * self._ds.metadata.unit["length"]["conv"] # Convert to code unit
        y = self.y * self._ds.metadata.unit["length"]["conv"] # Convert to code unit
        z = self.z * self._ds.metadata.unit["length"]["conv"] # Convert to code unit
        R = np.sqrt(x**2 + y**2 + z**2) / self._ds.metadata.unit["length"]["conv"] # Convert to user unit
        return R

    @property
    def p(self):
        r"""
        Particle absolute momentum.

        Notes
        -----
        p is calculated as follow :

        .. math::
            p = \sqrt{p_x^2+p_y^2+p_z^2}
        """
        px = self.px * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        py = self.py * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        pz = self.pz * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        p = np.sqrt(px**2 + py**2 + pz**2) / self._ds.metadata.unit["momentum"]["conv"] # Convert to user unit
        return p

    @property
    def thetax(self):
        r"""
        Particle polar angle for reference axis x (angle between px and the plane (py,pz)).

        Notes
        -----
        thetax is calculated as follow :

        .. math::
            \theta_x = \arccos{p_x/p}

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->py, y->pz, z->px).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system
        """
        px = self.px * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        p = self.p * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        thetax = np.arccos(px/p) / self._ds.metadata.unit["angle"]["conv"] # Convert to user unit
        return thetax

    @property
    def thetay(self):
        r"""
        Particle polar angle for reference axis y (angle between py and the plane (pz,px)).

        Notes
        -----
        thetay is calculated as follow :

        .. math::
            \theta_y = \arccos{p_y/p}

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->pz, y->px, z->py).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system
        """
        py = self.py * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        p = self.p * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        thetay = np.arccos(py/p) / self._ds.metadata.unit["angle"]["conv"] # Convert to user unit
        return  thetay

    @property
    def thetaz(self):
        r"""
        Particle polar angle for reference axis z (angle between px and the plane (py,pz)).

        Notes
        -----
        thetaz is calculated as follow :

        .. math::
            \theta = \arccos{p_x/p}

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->px, y->py, z->pz).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system
        """
        pz = self.pz * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        p = self.p * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        thetaz = np.arccos(pz/p) / self._ds.metadata.unit["angle"]["conv"] # Convert to user unit
        return thetaz

    @property
    def costhetax(self):
        r"""
        Particle polar angle for reference axis x (angle between px and the plane (py,pz)).

        Notes
        -----
        costhetax is calculated as follow :

        .. math::
            \cos(\theta_x) = p_x/p

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->py, y->pz, z->px).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system
        """
        px = self.px * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        p = self.p * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        costhetax = px/p # Dimensionless
        return costhetax

    @property
    def costhetay(self):
        r"""
        Particle polar angle for reference axis y (angle between py and the plane (pz,px)).

        Notes
        -----
        costhetay is calculated as follow :

        .. math::
            \cos(\theta_y) = p_y/p

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->pz, y->px, z->py).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system
        """
        py = self.py * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        p = self.p * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        costhetay = py/p # Dimensionless
        return costhetay

    @property
    def costhetaz(self):
        r"""
        Particle polar angle for reference axis z (angle between px and the plane (py,pz)).

        Notes
        -----
        costhetaz is calculated as follow :

        .. math::
            \cos(\theta_z) = p_z/p

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->px, y->py, z->pz).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system
        """
        pz = self.pz * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        p = self.p * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        costhetaz = pz/p # Dimensionless
        return costhetaz

    @property
    def phix(self):
        r"""
        Particle azimutal angle for reference axis x (angle between py and pz).

        Notes
        -----
        phix is calculated as follow :

        .. math::
            \phi_x = \arctan{p_z/p_y}

        with using the arctan2 function to have a result between -180 and 180 deg

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->py, y->pz, z->px).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system

        https://en.wikipedia.org/wiki/Atan2
        """
        pz = self.pz * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        py = self.py * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        phix = (np.arctan2(pz,py) + np.pi) / self._ds.metadata.unit["angle"]["conv"] # Convert to user unit
        return phix

    @property
    def phiy(self):
        r"""
        Particle azimutal angle for reference axis y (angle between pz and px).

        Notes
        -----
        phiy is calculated as follow :

        .. math::
            \phi_y = \arctan{p_x/p_z}

        with using the arctan2 function to have a result between -180 and 180 deg

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->pz, y->px, z->py).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system

        https://en.wikipedia.org/wiki/Atan2
        """
        px = self.px * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        pz = self.pz * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        phiy = (np.arctan2(px,pz) + np.pi) / self._ds.metadata.unit["angle"]["conv"] # Convert to user unit
        return phiy

    @property
    def phiz(self):
        r"""
        Particle azimutal angle for reference axis z (angle between px and py).

        Notes
        -----
        phiz is calculated as follow :

        .. math::
            \phi_z = \arctan{p_z/p_y}

        with using the arctan2 function to have a result between -180 and 180 deg

        References
        ----------
        https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations#To_spherical_coordinates
        with switching (x->px, y->py, z->pz).
        A figure can be found at https://en.wikipedia.org/wiki/Spherical_coordinate_system

        https://en.wikipedia.org/wiki/Atan2
        """
        py = self.py * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        px = self.px * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        phiz = (np.arctan2(py,px) + np.pi) / self._ds.metadata.unit["angle"]["conv"] # Convert to user unit
        return phiz

    @property
    def omegax(self):
        r"""
        Minimum solid angle in which the particle is contained (in sr).

        Notes
        -----
        omegax is calculated as follow :

        .. math::
            2 \pi (1 - \cos{\theta_x})

        References
        ----------
        https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
        """
        thetax = self.thetax * self._ds.metadata.unit["angle"]["conv"] # Convert to code unit
        omegax = 2*np.pi*(1-np.cos(thetax)) / self._ds.metadata.unit["solid_angle"]["conv"] # Convert to user unit
        return omegax

    @property
    def omegay(self):
        r"""
        Minimum solid angle in which the particle is contained (in sr).

        Notes
        -----
        omegay is calculated as follow :

        .. math::
            2 \pi (1 - \cos{\theta_y})

        References
        ----------
        https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
        """
        thetay = self.thetay * self._ds.metadata.unit["angle"]["conv"] # Convert to code unit
        omegay = 2*np.pi*(1-np.cos(thetay)) / self._ds.metadata.unit["solid_angle"]["conv"] # Convert to user unit
        return omegay

    @property
    def omegaz(self):
        r"""
        Minimum solid angle in which the particle is contained (in sr).

        Notes
        -----
        omegaz is calculated as follow :

        .. math::
            2 \pi (1 - \cos{\theta_z})

        References
        ----------
        https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
        """
        thetaz = self.thetaz * self._ds.metadata.unit["angle"]["conv"] # Convert to code unit
        omegaz = 2*np.pi*(1-np.cos(thetaz)) / self._ds.metadata.unit["solid_angle"]["conv"] # Convert to user unit
        return omegaz

    @property
    def etot(self):
        r"""
        Particle total energy.

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
        mass = self._ds.metadata.specie["mass"] # Already in code unit
        p = self.p * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        etot = np.sqrt(p**2 + mass**2) / self._ds.metadata.unit["energy"]["conv"] # Convert to user unit
        return etot

    @property
    def ekin(self):
        r"""
        Particle kinetic energy.

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
        mass = self._ds.metadata.specie["mass"] # Already in code unit
        etot = self.etot * self._ds.metadata.unit["energy"]["conv"] # Convert to code unit
        ekin = (etot - mass) / self._ds.metadata.unit["energy"]["conv"] # Convert to user unit
        return ekin

    @property
    def gamma(self):
        r"""
        Particle Lorentz factor.

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
        mass = self._ds.metadata.specie["mass"] # Already in code unit
        etot = self.etot * self._ds.metadata.unit["energy"]["conv"] # Convert to code unit
        gamma = self.etot/mass # Dimensionless
        return gamma

    @property
    def gammax(self):
        r"""
        Particle Lorentz factor along x axis.

        Notes
        -----
        gammax is calculated as follow :

        .. math::
            \gamma_x = \sqrt{1+(p_x/m_0)^2}

        with :math:`m_0` being the rest mass energy

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        mass = self._ds.metadata.specie["mass"] # Already in code unit
        px = self.px * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        gammax = np.sqrt(1+(self.px/mass)**2) # Dimensionless
        return gammax

    @property
    def gammay(self):
        r"""
        Particle Lorentz factor along y axis.

        Notes
        -----
        gammay is calculated as follow :

        .. math::
            \gamma_y = \sqrt{1+(p_y/m_0)^2}

        with :math:`m_0` being the rest mass energy

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        mass = self._ds.metadata.specie["mass"] # Already in code unit
        py = self.py * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        gammay = np.sqrt(1+(self.py/mass)**2) # Dimensionless
        return gammay

    @property
    def gammaz(self):
        r"""
        Particle Lorentz factor along z axis.

        Notes
        -----
        gammaz is calculated as follow :

        .. math::
            \gamma_z = \sqrt{1+(p_z/m_0)^2}

        with :math:`m_0` being the rest mass energy

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        mass = self._ds.metadata.specie["mass"] # Already in code unit
        pz = self.pz * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        gammaz = np.sqrt(1+(self.pz/mass)**2) # Dimensionless
        return gammaz

    @property
    def beta(self):
        r"""
        Particle normalized velocity.

        Notes
        -----
        beta is calculated as follow :

        .. math::
            \beta = \sqrt{1 - \frac{1}{\gamma^2}}

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        beta = np.sqrt(1.-1./self.gamma**2) # Dimensionless
        return beta

    @property
    def betax(self):
        r"""
        Particle normalized velocity along x axis.

        Notes
        -----
        betax is calculated as follow :

        .. math::
            \beta_x = \sqrt{1 - \frac{1}{\gamma_x^2}}

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        betax = np.sqrt(1.-1./self.gammax**2) # Dimensionless
        return betax

    @property
    def betay(self):
        r"""
        Particle normalized velocity along y axis.

        Notes
        -----
        betay is calculated as follow :

        .. math::
            \beta_y = \sqrt{1 - \frac{1}{\gamma_y^2}}

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        betay = np.sqrt(1.-1./self.gammay**2) # Dimensionless
        return betay

    @property
    def betaz(self):
        r"""
        Particle normalized velocity along z axis.

        Notes
        -----
        betaz is calculated as follow :

        .. math::
            \beta_z = \sqrt{1 - \frac{1}{\gamma_z^2}}

        References
        ----------
        https://en.wikipedia.org/wiki/Lorentz_factor
        """
        betaz = np.sqrt(1.-1./self.gammaz**2) # Dimensionless
        return betaz

    @property
    def ux(self):
        r"""
        Particle direction on x axis.

        Notes
        -----
        ux is calculated as follow :

        .. math::
            u_x = \frac{p_x}{p}
        """
        px = self.px * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        p = self.p * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        ux = px/p # Dimensionless
        return ux

    @property
    def uy(self):
        r"""
        Particle direction on y axis.

        Notes
        -----
        uy is calculated as follow :

        .. math::
            u_y = \frac{p_y}{p}
        """
        py = self.py * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        p = self.p * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        uy = py/p # Dimensionless
        return uy

    @property
    def uz(self):
        r"""
        Particle direction on z axis.

        Notes
        -----
        uz is calculated as follow :

        .. math::
            u_z = \frac{p_z}{p}
        """
        pz = self.pz * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        p = self.p * self._ds.metadata.unit["momentum"]["conv"] # Convert to code unit
        uz = pz/p # Dimensionless
        return uz

    @property
    def m(self):
        r"""
        Particle relativistic mass.

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
        mass = self._ds.metadata.specie["mass"] # Already in code unit
        gamma = self.gamma # Dimensionless
        m = gamma * m / self._ds.metadata.unit["momentum"]["conv"] # Convert to user unit
        return m

    @property
    def wekin(self):
        r"""
        Particle kinetic energy density.

        Notes
        -----
        wekin is calculated as follow :

        .. math::
            E_{kin} \times w
        """
        wekin = self.w * self.ekin
        return wekin

    # @property
    # def womegax(self): ???
    #     r"""
    #     Particle kinetic energy density.
    #
    #     Notes
    #     -----
    #     wekin is calculated as follow :
    #
    #     .. math::
    #         E_{kin} \times w
    #     """
    #     return self.w * self.ekin
