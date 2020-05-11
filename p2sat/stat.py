#coding:utf8

"""
Get global statistics from datasets.
"""

import numpy as _np
from .datasets import PhaseSpace as _PhaseSpace

def expected_value(ds,qty,select=None):
    r"""
    Returns expected value of given quantity.

    Parameters
    ----------
    ds : {PhaseSpace, ScalarField, EventLocation}
        Dataset to use.
    qty : str or np.array
        Quantity to consider.
    select : dict, optional
        Filtering dictionary.

    Notes
    -----
    expected_value is defined as sum(p*qty) with p=w/sum(w)
    """
    if isinstance(ds, _PhaseSpace):
        r = ds.read
        qty = r.quantity(qty, select=select)
        w = r.quantity('w',select=select)

        p = w/sum(w)

        return sum(p*qty)
    else:
        raise NotImplementedError()

def variance(ds,qty,select=None):
    r"""
    Returns variance of given quantity.

    Parameters
    ----------
    ds : {PhaseSpace, ScalarField, EventLocation}
        Dataset to use.
    qty : str or np.array
        Quantity to consider.
    select : dict, optional
        Filtering dictionary.

    Notes
    -----
    variance is defined as expected_value((qty - expected_value(qty))**2)
    """
    if isinstance(ds, _PhaseSpace):
        r = ds.read
        qty = r.quantity(qty,select=select)
        w = r.quantity("w",select=select)

        ev = ds.expected_value

        var = ev((qty - ev(qty))**2)

        return var
    else:
        raise NotImplementedError()

def standard_deviation(ds,qty,select=None):
    r"""
    Returns standard deviation of given quantity.

    Parameters
    ----------
    ds : {PhaseSpace, ScalarField, EventLocation}
        Dataset to use.
    qty : str or np.array
        Quantity to consider.
    select : dict, optional
        Filtering dictionary.

    Notes
    -----
    standard_deviation is defined as the square root of variance
    """
    return _np.sqrt(variance(ds, qty, select=select))

def covariance(ds,qty1,qty2,select=None):
    r"""
    Returns covariance of given quantities.

    Parameters
    ----------
    ds : {PhaseSpace, ScalarField, EventLocation}
        Dataset to use.
    qty1 : str or np.array
        Quantity to consider.
    qty2 : str or np.array
        Quantity to consider
    select : dict, optional
        Filtering dictionary.

    Notes
    -----
    covariance is defined as expected_value((qty1-expected_value(qty1)) * (qty2-expected_value(qty2)))
    """
    if isinstance(ds, _PhaseSpace):
        r = ds.read
        qty1 = r.quantity(qty1,select=select)
        qty2 = r.quantity(qty2,select=select)

        ev = ds.expected_value

        cov = ev((qty1-ev(qty1))*(qty2-ev(qty2)))

        return cov
    else:
        raise NotImplementedError()

def correlation_coefficient(ds,qty1,qty2,select=None):
    r"""
    Returns correlation coefficient of given quantities.

    Parameters
    ----------
    ds : {PhaseSpace, ScalarField, EventLocation}
        Dataset to use.
    qty1 : str or np.array
        Quantity to consider.
    qty2 : str or np.array
        Quantity to consider
    select : dict, optional
        Filtering dictionary.

    Notes
    -----
    correlation_coefficient is defined as covariance(qty1,qty2)/(standard_deviation(qty1)*standard_deviation(qty2))
    """
    std = ds.standard_deviation

    cc = covariance(ds, qty1, qty2, select=select)/(std(ds, qty1, select=select)*std(ds, qty2, select=select))

    return cc

def total_number(ds, select=None):
    r"""
    Return total number of particles contained in the PhaseSpace.

    Parameters
    ----------
    ds : PhaseSpace
        Dataset to use.
    select : dict, optional
        Filtering dictionary.

    See Also
    --------
    data.select
    """
    if isinstance(ds, _PhaseSpace):
        return sum(ds.read.quantity("w", select=select))
    else:
        raise NotImplementedError()

def total_energy(ds, unit="J", select=None):
    r"""
    Return total energy contained in the dataset.

    Parameters
    ----------
    ds : {PhaseSpace, ScalarField}
        Dataset to use.
    unit : str, optional
        unit of energy. Available. are 'J' and 'MeV'. Default is 'J'
    select : dict, optional
        Filtering dictionary.

    See Also
    --------
    data.select
    """
    if isinstance(ds, _PhaseSpace):
        r = ds.read
        w = r.quantity('w',select=select)
        ekin = r.quantity('ekin',select=select)

        E_MeV = 1e6 * sum(w*ekin)/ds.metadata.unit["energy"]["conv"]
        E_J = E_MeV * 1e6 * 1.6e-19

        if unit == "MeV":
            return E_MeV
        elif unit =="J":
            return E_J
        else:
            raise NameError("Unknown unit name.")
    else:
        raise NotImplementedError()

def total_charge(ds, unit="C", select=None):
    r"""
    Return total charge of the PhaseSpace.

    Parameters
    ----------
    ds : PhaseSpace
        Dataset to use.
    unit : str, optional
        unit of charge. Available are 'C' and 'nC', 'pC'. Default is 'nC'.
    select : dict, optional
        Filtering dictionary.

    See Also
    --------
    data.select
    """
    if isinstance(ds, _PhaseSpace):
        Q = ds.metadata.specie["charge"] * 1.6e-19 * sum(ds.read.quantity("w", select=select))
        if unit == "C":
            return Q
        elif unit =="nC":
            return 1e9 * Q
        elif unit =="pC":
            return 1e12 * Q
        else:
            raise NameError("Unknown unit name.")
    else:
        raise NotImplementedError()

def peak_brilliance(ds):
    r"""
    """
    raise NotImplementedError()
