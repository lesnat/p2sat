#coding:utf8
import numpy as np

class _Stat(object):
    """
    Get global statistics from phase space data.
    """
    def __init__(self,PhaseSpace):
        self._ps = PhaseSpace

    def expected_value(self,axis,select=None):
        """
        Returns expected value of given axis.

        Parameters
        ----------
        axis : str or np.array
            axis to consider
        select : dict, optional
            filtering dictionary

        Returns
        -------
        E : float
            expected value of given axis

        Notes
        -----
        expected_value is defined as sum(p*axis) with p=w/sum(w)
        """
        d=self._ps.data
        axis = d.get_axis(axis,select=select)
        w = d.get_axis('w',select=select)

        p = w/sum(w)

        return sum(p*axis)

    def variance(self,axis,select=None):
        """
        Returns variance of given axis.

        Parameters
        ----------
        axis : str or np.array
            axis to consider
        select : dict, optional
            filtering dictionary

        Returns
        -------
        var : float
            variance of given axis

        Notes
        -----
        variance is defined as expected_value((axis - expected_value(axis))**2)
        """
        d=self._ps.data
        axis = d.get_axis(axis,select=select)
        w = d.get_axis("w",select=select)

        ev = self.expected_value

        var = ev((axis - ev(axis))**2)

        return var

    def standard_deviation(self,axis,select=None):
        """
        Returns standard deviation of given axis.

        Parameters
        ----------
        axis : str or np.array
            axis to consider
        select : dict, optional
            filtering dictionary

        Returns
        -------
        std : float
            standard deviation of given axis

        Notes
        -----
        standard_deviation is defined as the square root of variance
        """
        return np.sqrt(self.variance(axis,select=select))

    def covariance(self,axis1,axis2,select=None):
        """
        Returns covariance of given axes.

        Parameters
        ----------
        axis1 : str or np.array
            axis to consider
        axis2 : str or np.array
            axis to consider
        select : dict, optional
            filtering dictionary

        Returns
        -------
        cov : float
            covariance of given axes

        Notes
        -----
        covariance is defined as expected_value((axis1-expected_value(axis1)) * (axis2-expected_value(axis2)))
        """
        d=self._ps.data
        axis1 = d.get_axis(axis1,select=select)
        axis2 = d.get_axis(axis2,select=select)

        ev = self.expected_value

        cov = ev((axis1-ev(axis1))*(axis2-ev(axis2)))

        return cov

    def correlation_coefficient(self,axis1,axis2,select=None):
        """
        Returns correlation coefficient of given axes.

        Parameters
        ----------
        axis1 : str or np.array
            axis to consider
        axis2 : str or np.array
            axis to consider
        select : dict, optional
            filtering dictionary

        Returns
        -------
        cc : float
            correlation coefficient of given axes

        Notes
        -----
        correlation_coefficient is defined as covariance(axis1,axis2)/(standard_deviation(axis1)*standard_deviation(axis2))
        """
        std = self.standard_deviation

        cc = self.covariance(axis1,axis2,select=select)/(std(axis1,select=select)*std(axis2,select=select))

        return cc

    def total_energy(self,unit="J",select=None):
        """
        Return total energy contained in the phase space

        Parameters
        ----------
        unit : str, optional
            unit of energy. Available are 'J' and 'MeV'. Default is 'J'
        select : dict, optional
            filtering dictionary

        See Also
        --------
        data.select
        """
        d = self._ps.data
        w = d.get_axis('w',select=select)
        ekin = d.get_axis('ekin',select=select)

        E_MeV = sum(w*ekin)
        E_J = E_MeV * 1e6 * 1.6e-19

        if unit == "MeV":
            return E_MeV
        elif unit =="J":
            return E_J
        else:
            raise NameError("Unknown unit name.")
