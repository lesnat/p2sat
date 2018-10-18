#coding:utf8
import numpy as np

class _Hist(object):
    """
    Create histograms from raw data.
    """
    def __init__(self,PhaseSpace):
        self._ps=PhaseSpace

    def hn(self,axis,weight="w",bwidth=None,brange=None,normed=True,select=None):
        """
        Create and return the n-dimensional histo of axis list.

        Parameters
        ----------
        axis : list of str/np.array
            list of axis to hist
        bwidth : list of float, optional
            list of bin width. If a bwidth element is None, a calculation is done to have 100 bins in the correspondant axis
        brange : list of list of 2 float, optional
            list of bin minimum and maximum. If a brange element is None, the minimum/maximum of the axis is taken
        normed : bool or list of bool, optional
            weight normalization. If a normed element is True, the bin width is taken as weight normalization
        select : dict, optional
            filtering dictionary

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
        and value/range of filtering axis as dict values.

        Examples
        --------
        Using `eps` as a `PhaseSpace` instance

        >>> w,x = eps.hist.hn(['x'],
        ...                   bwidth=[50],brange=[[0,1000]],
        ...                   normed=[True],
        ...                   select={'ekin':(0.511,None)})
        ...

        returns the number of particles with :math:`ekin \in [0.511, +\infty] MeV` in function of x
        normed=[True] to divide weight by bin width, so weight unit is Number/um

        >>> w,r,ekin=eps.hist.hn(['r','ekin'],
        ...                      bwidth=[10.0,0.1],
        ...                      brange=[[0,1000],[0.1,50.0]],
        ...                      select={'x':150})
        ...

        returns the number of particle per um per MeV at x=150 um

        See Also
        --------
        data.select
        """
        # Get a shortcut to data object
        d=self._ps.data

        # # Get a copy of the axes
        for i,ax in enumerate(axis):
            axis[i] = d.get_axis(ax,select=select)

        # Get weight array
        w   = d.get_axis(weight,select=select)

        # Define bin range
        if brange is None   : brange=[[None,None]]*len(axis)
        for i,_ in enumerate(axis):
            if brange[i][0] is None:brange[i][0]=min(axis[i])
            if brange[i][1] is None:brange[i][1]=max(axis[i])

        # Define bin width
        if bwidth is None   : bwidth=[None]*len(axis)
        blen=[None]*len(axis)
        for i,_ in enumerate(axis):
            if bwidth[i] is not None:
                # Bin width is over-estimated to fit with bin range
                blen[i]   = int(np.ceil((brange[i][1] + bwidth[i] - brange[i][0])/bwidth[i]))
            else:
                # Default is 100 bins
                blen[i]   = 100
                bwidth[i] = float(brange[i][1] - brange[i][0])/blen[i]

        # Calculate bins
        bins=[]
        for i,_ in enumerate(axis):
            bins.append(np.linspace(brange[i][0],brange[i][1],blen[i]))

        # Get weight normalization
        if type(normed) is bool:
            if normed:
                normed=[True]*len(axis)
            else:
                normed=[False]*len(axis)
        wnorm = 1.
        for i,_ in enumerate(axis):
            if normed[i]: wnorm*=bwidth[i]

        # Calculate the multi dimensional histo, normalized by wnorm
        h,b=np.histogramdd(axis,weights=w/wnorm,bins=bins)

        # Return the bins and histo
        return b,h

    def h1(self,axis,weight="w",bwidth=None,brange=None,normed=True,select=None):
        """
        Create and return the 1 dimensional histogram of given axis.

        Parameters
        ----------
        axis : str or np.array
            axis to hist
        bwidth : float, optional
            bin width. If None, a calculation is done to have 10 bins in the axis
        brange : list of 2 float, optional
            bin maximum and minimum. If a brange element is None, the axis minimum/maximum is taken
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
        the h1 method is just a different way to call the generic method hn

        See Also
        --------
        hist.hn, hist.h2, hist.h3
        """
        if not brange : brange = [None,None]

        b,h=self.hn([axis],weight=weight,bwidth=[bwidth],brange=[brange],normed=normed,select=select)

        return b[0],h

    def h2(self,axis1,axis2,weight="w",bwidth1=None,bwidth2=None,brange1=None,brange2=None,normed=True,select=None):
        """
        Create and return the 2 dimensional histogram of given axis.

        Parameters
        ----------
        axis1,axis2 : str or np.array
            axis to hist
        bwidth1,bwidth2 : float, optional
            bin width. If None, a calculation is done to have 10 bins in the axis
        brange1,brange2 : list of 2 float, optional
            bin maximum and minimum. If a brange element is None, the axis minimum/maximum is taken
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
        the h2 method is just a different way to call the generic method hn

        See Also
        --------
        hist.hn, hist.h1, hist.h3
        """
        if not brange1 : brange1 = [None,None]
        if not brange2 : brange2 = [None,None]

        b,h=self.hn([axis1,axis2],weight=weight,bwidth=[bwidth1,bwidth2],brange=[brange1,brange2],normed=normed,select=select)

        return b[0],b[1],h

    def h3(self,axis1,axis2,axis3,weight="w",bwidth1=None,bwidth2=None,bwidth3=None,brange1=None,brange2=None,brange3=None,normed=True,select=None):
        """
        Create and return the 3 dimensional histogram of given axis.

        Parameters
        ----------
        axis1,axis2,axis3 : str or np.array
            axis to hist
        bwidth1,bwidth2,bwidth3 : float, optional
            bin width. If None, a calculation is done to have 10 bins in the axis
        brange1,brange2,brange3 : list of 2 float, optional
            bin maximum and minimum. If a brange element is None, the axis minimum/maximum is taken
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
        the h3 method is just a different way to call the generic method hn

        See Also
        --------
        hist.hn, hist.h1, hist.h2
        """
        if not brange1 : brange1 = [None,None]
        if not brange2 : brange2 = [None,None]
        if not brange3 : brange3 = [None,None]

        b,h=self.hn([axis1,axis2,axis3],weight="w",bwidth=[bwidth1,bwidth2,bwidth3],brange=[brange1,brange2,brange3],normed=normed,select=select)

        return b[0],b[1],b[2],h

    def f1(self,axis,func_name,weight="w",return_fit=False,verbose=True,**kargs):
        """
        Fit a 1D histogram with given law.

        Parameters
        ----------
        axis : str or np.array
            axis to fit
        func_name : str
            name of the fit law. Available are `exp` for exponential law and `gauss` for gaussian law
        return_fit : bool, optional
            returns the spectrum instead of fited parameters
        verbose : bool, optional
            verbosity
        kargs : dict, optional
            dictionnary to pass to the hist.h1 method

        Returns
        -------
        x : np.array
            fit abscissa
        param1,param2 : float,optional
            fit parameters. Returned if `return_fit=False` (default)
        w : np.array, optional
            fit weight. Returned if `return_fit=True`

        Notes
        -----
        The `exp` law is defined as
        :math:`\\frac{N}{T} \exp{(-x/T)}`
        and returns fit parameters N,T.

        The `gauss` law is defined as
        :math:`\\frac{N}{ \sigma \sqrt{2 \pi}} \exp{(-\\frac{(x-\mu)^2}{2 \sigma^2})}`
        and returns fit parameters N,sigma,mu.
        """
        # Get the hist data
        if verbose:print("Processing 1D \"{}\" fit of \"{}\" ...".format(func_name,str(axis)))
        x,w = self._ps.hist.h1(axis,**kargs)

        bwidth = np.array([x[i+1]-x[i] for i in range(len(w))])
        Ntot = sum(w*bwidth)
        # Define fit function and default values for fit parameters
        if func_name=="exp":
            f = lambda x,N,T: N*np.exp(-x/T)/T
            p0 = [Ntot,1]
            param = ["N","T"]
        elif func_name=="gauss":
            f = lambda x,N,sigma,mu: N/(np.sqrt(2*np.pi) * sigma) * np.exp(-(x-mu)**2/(2*sigma**2))
            p0 = [Ntot,x.std(),0]
            param = ["N","sigma","mu"]
        else:
            raise NameError("Unknown func_name.")

        # Fit the curve
        from scipy.optimize import curve_fit
        popt,pcov = curve_fit(f,x[:-1],w,p0=p0)
        perr = np.sqrt(np.diag(pcov)) # Estimated error
        perr_pc = (1. - (popt - perr)/popt) * 100

        # Print estimated errors
        if verbose:
            print('Fit parameters :')
            for i,p in enumerate(param):
                try:
                    if p=="N": raise TypeError
                    unit = self._ps.data.raw.units[axis]
                except TypeError:
                    unit = ""
                print(u"    {} = {: .4E} ± {:.4E} {:s} ({:.2F} %)".format(p.ljust(7),popt[i],perr[i], unit.ljust(5), perr_pc[i]))

        # Choice of the return values
        if return_fit:
            # Return axis and spectrum
            return x,f(x,*popt)
        else:
            # Format the result in a list
            res = [x]
            for e in popt:
                res.append(e)
            # Return axis and fit parameters
            return res
