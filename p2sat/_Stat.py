#coding:utf8
import numpy as np

class _Stat(object):
  """
  Allows to do statistics with p2sat data
  """
  def __init__(self,PhaseSpace):
    self._ps = PhaseSpace

  def expected_value(self,axis,p=None,select=None):
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
    with w being the statistical weight of the configuration

    """
    r=self._ps.data
    if type(axis) is str:axis=eval("r.%s"%axis)

    w = np.array(r.w)
    if select is not None:
      w = r.select(w,faxis=select.keys(),frange=select.values())
      axis = r.select(axis,faxis=select.keys(),frange=select.values())

    if p is None:
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
    r=self._ps.data
    if type(axis) is str:axis=eval("r.%s"%axis)

    p=None
    if select is not None:
      axis = r.select(axis,faxis=select.keys(),frange=select.values())
      w = r.select(r.w,faxis=select.keys(),frange=select.values())
      p = w/sum(w)

    ev = self.expected_value

    var = ev((axis - ev(axis,p=p))**2,p=p)

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
    r=self._ps.data
    if type(axis1) is str:axis1=eval("r.%s"%axis1)
    if type(axis2) is str:axis2=eval("r.%s"%axis2)

    p = None
    if select is not None:
      axis1 = r.select(axis1,faxis=select.keys(),frange=select.values())
      axis2 = r.select(axis2,faxis=select.keys(),frange=select.values())
      w = r.select(r.w,faxis=select.keys(),frange=select.values())
      p = w/sum(w)

    ev = self.expected_value

    cov = ev((axis1-ev(axis1,p=p))*(axis2-ev(axis2,p=p)),p=p)

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
