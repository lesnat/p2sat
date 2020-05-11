#coding:utf8
"""
Create histograms from raw data.
"""

import numpy as _np

def histNd(ds, qty, weight="w", bwidth=None, brange=None, normed=True, select=None):
    r"""
    Create and return the n-dimensional histo of qty list.

    Parameters
    ----------
    ds : {PhaseSpace, EventLocation}
        Dataset to use.
    qty : list of str/np.array
        List of qty to hist.
    bwidth : list of float, optional
        List of bin width. If a bwidth element is None, a calculation is done to have 100 bins in the correspondant qty.
    brange : list of list of 2 float, optional
        List of bin minimum and maximum. If a brange element is None, the minimum/maximum of the qty is taken.
    normed : bool or list of bool, optional
        Weight normalization. If a normed element is True, the bin width is taken as weight normalization.
    select : dict, optional
        Filtering dictionary.

    Returns
    -------
    b : np.array
        bins
    h : np.array
        number of particles per bin unit

    Notes
    -----
    The weight normalization allows to be independant of choosen bin width.
    If normalized, the weight unit is `Number/(unit1 x unit2 x ...)` with
    `unit1, unit2, ...` the units of axes `1, 2, ...`.

    If the given maximum bin range does not match with an int number of bins, the bin width is over sized.

    The select dictionary takes the name of filtering axes as dict keys,
    and value/range of filtering qty as dict values.

    Examples
    --------
    >>> eps = ExamplePhaseSpace()
    >>> w,x = eps.hist.histNd(['x'],
    ...                   bwidth=[50],brange=[[0,1000]],
    ...                   normed=[True],
    ...                   select={'ekin':(0.511,None)})
    ...

    returns the number of particles with :math:`ekin \in [0.511, +\infty] MeV` in function of x
    normed=[True] to divide weight by bin width, so weight unit is Number/um

    >>> w,r,ekin=eps.hist.histNd(['r','ekin'],
    ...                      bwidth=[10.0,0.1],
    ...                      brange=[[0,1000],[0.1,50.0]],
    ...                      select={'x':150})
    ...

    returns the number of particle per um per MeV at x=150 um

    See Also
    --------
    data.select
    """
    # Get a copy of the axes
    for i, q in enumerate(qty):
        qty[i] = ds.read.quantity(q, select=select)

    # Get weight array
    w   = ds.read.quantity(weight,select=select)

    # Default bin range and width
    if brange is None   : brange=[[None,None]]*len(qty)
    if bwidth is None   : bwidth=[None]*len(qty)
    if type(normed) is bool: normed=[normed] * len(qty)

    # Define weight normalization
    wnorm = 1.

    # Calculate bins
    bins=[]
    for i,_ in enumerate(qty):
        # Default bin range are min and max values
        if brange[i][0] is None: brange[i][0] = min(qty[i])
        if brange[i][1] is None: brange[i][1] = max(qty[i])
        # If min == max, change brange to make an histogram anyway
        if brange[i][0] == brange[i][1]:
            # If the values are 0, the last bin will be 0.1. It will be 2 times the qty value otherwise.
            if brange[i][0] == 0.:
                brange[i][1] = 0.1
            else:
                brange[i][0] = 0.
                brange[i][1] *= 2.

        # Default bin width
        if bwidth[i] is None: bwidth[i] = (brange[i][1] - brange[i][0])/100

        # Calculate bin edges
        b = _np.arange(brange[i][0], brange[i][1]+2*bwidth[i], bwidth[i]) # max is excluded
        bins.append(b)

        # Normalize weight/bin
        if normed[i]: wnorm*=bwidth[i]

    # Calculate the multi dimensional histo, normalized by wnorm
    h,b=_np.histogramdd(qty,weights=w/wnorm,bins=bins)

    # Return the bins and histo
    return b,h

def hist1d(ds, qty, weight="w", bwidth=None, brange=None, normed=True, select=None):
    r"""
    Create and return the 1 dimensional histogram of given quantity.

    Parameters
    ----------
    ds : {PhaseSpace, EventLocation}
        Dataset to use.
    qty : str or np.array
        qty to hist
    bwidth : float, optional
        bin width. If None, a calculation is done to have 10 bins in the qty
    brange : list of 2 float, optional
        bin maximum and minimum. If a brange element is None, the qty minimum/maximum is taken
    normed : bool, optional
        weight normalization. If a normed element is True, the bin width is taken as weight normalization
    select : dict, optional
        filtering dictionnary

    Returns
    -------
    b : np.array
        bins
    h : np.array
        histogram

    Notes
    -----
    the hist1d method is just a different way to call the generic method histNd

    See Also
    --------
    hist.histNd, hist.hist2d, hist.hist3d
    """
    if not brange : brange = [None,None]

    b,h=histNd(ds,[qty],weight=weight,bwidth=[bwidth],brange=[brange],normed=normed,select=select)

    return b[0],h

def hist2d(ds, qty1, qty2, weight="w", bwidth1=None, bwidth2=None, brange1=None, brange2=None, normed=True, select=None):
    r"""
    Create and return the 2 dimensional histogram of given quantity.

    Parameters
    ----------
    ds : {PhaseSpace, EventLocation}
        Dataset to use.
    qty1,qty2 : str or np.array
        qty to hist
    bwidth1,bwidth2 : float, optional
        bin width. If None, a calculation is done to have 10 bins in the qty
    brange1,brange2 : list of 2 float, optional
        bin maximum and minimum. If a brange element is None, the qty minimum/maximum is taken
    normed : bool, optional
        weight normalization. If a normed element is True, the bin width is taken as weight normalization
    select : dict, optional
        filtering dictionnary

    Returns
    -------
    b1,b2 : np.array
        bins
    h : np.array
        histogram

    Notes
    -----
    the hist2d method is just a different way to call the generic method histNd

    See Also
    --------
    hist.histNd, hist.hist1d, hist.hist3d
    """
    if not brange1 : brange1 = [None,None]
    if not brange2 : brange2 = [None,None]

    b,h=histNd(ds,[qty1,qty2],weight=weight,bwidth=[bwidth1,bwidth2],brange=[brange1,brange2],normed=normed,select=select)

    return b[0],b[1],h

def hist3d(ds, qty1, qty2, qty3, weight="w", bwidth1=None, bwidth2=None, bwidth3=None, brange1=None, brange2=None, brange3=None, normed=True, select=None):
    r"""
    Create and return the 3 dimensional histogram of given quantity.

    Parameters
    ----------
    ds : {PhaseSpace, EventLocation}
        Dataset to use.
    qty1,qty2,qty3 : str or np.array
        qty to hist
    bwidth1,bwidth2,bwidth3 : float, optional
        bin width. If None, a calculation is done to have 10 bins in the qty
    brange1,brange2,brange3 : list of 2 float, optional
        bin maximum and minimum. If a brange element is None, the qty minimum/maximum is taken
    normed : bool, optional
        weight normalization. If a normed element is True, the bin width is taken as weight normalization
    select : dict, optional
        filtering dictionnary

    Returns
    -------
    b1,b2,b3 : np.array
        bins
    h : np.array
        histogram

    Notes
    -----
    the hist3d method is just a different way to call the generic method histNd

    See Also
    --------
    hist.histNd, hist.hist1d, hist.hist2d
    """
    if not brange1 : brange1 = [None,None]
    if not brange2 : brange2 = [None,None]
    if not brange3 : brange3 = [None,None]

    b,h=histNd(ds,[qty1,qty2,qty3],weight="w",bwidth=[bwidth1,bwidth2,bwidth3],brange=[brange1,brange2,brange3],normed=normed,select=select)

    return b[0],b[1],b[2],h

def fit1d(ds, qty, f, weight="w", p0=None, verbose=True, **kargs):
    r"""
    Fit a 1D histogram with given function.

    Parameters
    ----------
    ds : PhaseSpace
        Dataset to use.
    qty : str or np.array
        Quantity to fit.
    f : function
        Fitting function.
    verbose : bool, optional
        Verbosity
    kargs : dict, optional
        Dictionnary to pass to the hist.hist1d method.

    Returns
    -------
    popt : float
        fit parameters

    See Also
    --------
    scipy.optimize.curve_fit
    """
    if verbose:pass
    # Get the hist data
    x, w = hist1d(ds, qty, weight=weight, **kargs)

    # Fit the curve
    from scipy.optimize import curve_fit
    popt, pcov  = curve_fit(f, x[:-1], w, p0=p0)
    perr        = _np.sqrt(_np.diag(pcov)) # Estimated error
    perr_pc     = (1. - (popt - perr)/popt) * 100

    # Print estimated errors
    if verbose:
        pass

    return popt
