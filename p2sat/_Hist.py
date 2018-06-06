#coding:utf8
import numpy as np

class _Hist(object):
  """
  Histograms
  """
  def __init__(self,PhaseSpace):
    self._ps=PhaseSpace
    r=self._ps.raw

  # def hn(self,axis,blen=None,bwidth=None,brange=None,wnorm=None,select=None):
  #   """
  #   Create and return the n-dimensional histo of axis list.
  #
  #   Parameters
  #   ----------
  #   axis : list of str/np.array
  #     list of axis to hist
  #   blen : list of int, optional
  #     list of number of bins. If a blen element is None, the default value is 10
  #   bwidth : list of float, optional
  #     list of bin width. If a bwidth element is None, a calculation is done to have 10 bins in the correspondant axis
  #   brange : list of list of 2 float, optional
  #     list of bin minimum and maximum. If a brange element is None, the minimum/maximum of the axis is taken
  #   wnorm : list of float, optional
  #     weight normalization. If a wnorm element is None, the bin width is taken
  #   select : dict, optional
  #     filtering dictionary
  #
  #   Returns
  #   -------
  #   b : np.array
  #     bins
  #   h : np.array
  #     number of particles per bin unit
  #
  #   Notes
  #   -----
  #   TODO: If the given maximum bin range does not match with an int number of bins, the last bin is oversized ??
  #   it reduce bwidth to fit brange[1]-brange[0] with a int nb of bins
  #
  #   If blen and bwidth are both defined, priority is given to blen.
  #
  #   Examples
  #   --------
  #
  #   >>> hn(['x'],bwidth=[50],brange=[[0,1000]],wnorm=[1.0],select={'ekin':(0.511,None)})
  #
  #   returns the number of particles with :math:`ekin \in [0.511, +\infty] MeV` in function of x
  #   wnorm=[1.0] to not divide nb of particles by bin width (otherwise number per um)
  #
  #   >>> hn(['r','ekin'],bwidth=[10.0,0.1],brange=[[0,1000],[0.1,50.0]],select={'x':150})
  #
  #   returns a number of e- per um per MeV at x=150 um
  #   """
  #   r=self._ps.raw
  #
  #   # Get a copy the axis from a str if needed
  #   for i,ax in enumerate(axis):
  #     if type(ax) is str:ax = eval("r.%s"%ax)
  #     axis[i]=np.array(ax)
  #
  #   # Get a copy of particle statistical weight
  #   w   = np.array(r.w)
  #
  #   # Filter the data if needed
  #   if select is not None:
  #     for faxis,frange in select.items():
  #       w = r.select(w,faxis,frange) # FIXME: bug if more than 1 axis
  #       for i,ax in enumerate(axis):
  #         axis[i]=r.select(ax,faxis,frange)
  #
  #   # Define default bin range
  #   if not brange: brange=[[None,None]]*len(axis)
  #   for i,br in enumerate(brange):
  #     if br[0] is None:brange[i][0]=min(axis[i])
  #     if br[1] is None:brange[i][1]=max(axis[i])
  #   """
  #   # Define default bin width
  #
  #   # Calculate blen from bwidth
  #
  #   # Define default number of bins
  #   if not blen: blen=[None]*len(axis)
  #   for i,bl in enumerate(blen):
  #     if bl is None: blen[i]=int(np.ceil((brange[i][1] + bwidth[i] - brange[i][0])/bwidth[i]))
  #
  #
  #   bins.append(np.linspace(brange[i][0],brange[i][1]+bwidth[i],blen[i]) # FIXME: nope because of bwidth possibly not def.
  #   """
  #   # Calculate number of bins from bwidth if needed
  #   if not bwidth: bwidth=[None]*len(axis)
  #   for i,bw in enumerate(bwidth):
  #     blen=brange[i][1] - brange[i][0]
  #     if bw is None: bwidth[i]=blen/10.
  #
  #   # Define default bin width + number of bins
  #   nbins=[]
  #   if not bwidth: bwidth=[None]*len(axis)
  #   for i,bw in enumerate(bwidth):
  #     blen=brange[i][1] - brange[i][0]
  #     if bw is None: bwidth[i]=blen/10.
  #     nbins.append(int(np.ceil(np.nan_to_num(blen/bwidth[i])))) # TODO: delete nan_to_num?
  #
  #   # Construct the bins list
  #   bins=[]
  #   for i,ax in enumerate(axis):
  #     if nbins[i]>0:
  #       bins.append(np.linspace(brange[i][0],brange[i][1]+bwidth[i],nbins[i]+1))
  #     else:
  #       bins.append(1)
  #
  #     #bins.append(np.arange(brange[i][0],brange[i][1],bwidth[i]))
  #
  #   # Define weight normalization
  #   if not wnorm: wnorm=[None]*len(axis)
  #   for i,wn in enumerate(wnorm):
  #     if wn is None: wnorm[i]=bwidth[i]
  #
  #   # Calculate the multi dimensional histo, normalized by wnorm
  #   h,b=np.histogramdd(axis,weights=w/np.product(wnorm),bins=bins)
  #
  #   # Return the bins and histo
  #   return b,h

  def hn(self,axis,blen=None,bwidth=None,brange=None,wnorm=None,select=None):
    """
    Create and return the n-dimensional histo of axis list.

    Parameters
    ----------
    axis : list of str/np.array
      list of axis to hist
    blen : list of int, optional
      list of number of bins. If a blen element is None, the default value is 10
    bwidth : list of float, optional
      list of bin width. If a bwidth element is None, a calculation is done to have 10 bins in the correspondant axis
    brange : list of list of 2 float, optional
      list of bin minimum and maximum. If a brange element is None, the minimum/maximum of the axis is taken
    wnorm : list of float, optional
      weight normalization. If a wnorm element is None, the bin width is taken
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
    TODO: If the given maximum bin range does not match with an int number of bins, the last bin is oversized ??
    it reduce bwidth to fit brange[1]-brange[0] with a int nb of bins

    If blen and bwidth are both defined, priority is given to blen.

    Examples
    --------

    >>> hn(['x'],bwidth=[50],brange=[[0,1000]],wnorm=[1.0],select={'ekin':(0.511,None)})

    returns the number of particles with :math:`ekin \in [0.511, +\infty] MeV` in function of x
    wnorm=[1.0] to not divide nb of particles by bin width (otherwise number per um)

    >>> hn(['r','ekin'],bwidth=[10.0,0.1],brange=[[0,1000],[0.1,50.0]],select={'x':150})

    returns a number of e- per um per MeV at x=150 um
    """
    r=self._ps.raw

    # Get a copy the axis from a str if needed
    for i,ax in enumerate(axis):
      if type(ax) is str:ax = eval("r.%s"%ax)
      axis[i]=np.array(ax)

    # Get a copy of particle statistical weight
    w   = np.array(r.w)

    # Filter the data if needed
    # if select is not None:
    #   for faxis,frange in select.items():
    #     print(len(w))
    #     w = r.select(w,faxis,frange) # FIXME: bug if more than 1 axis
    #     for i,ax in enumerate(axis):
    #       print(len(ax))
    #       axis[i]=r.select(ax,faxis,frange)
    if select is not None:
      w = r.select(w,faxis=select.keys(),frange=select.values())
      for i,ax in enumerate(axis):
        axis[i]=r.select(ax,faxis=select.keys(),frange=select.values())

    # Define default bin range
    if wnorm is None    : wnorm=[None]*len(axis)
    if brange is None   : brange=[[None,None]]*len(axis)
    if bwidth is None   : bwidth=[None]*len(axis)
    if blen is None     : blen=[None]*len(axis)
    bins=[]
    for i,_ in enumerate(axis):
      if brange[i][0] is None:brange[i][0]=min(axis[i])
      if brange[i][1] is None:brange[i][1]=max(axis[i])
    #   if brange[i][0]==brange[i][1]:
    #     brange[i][1]=2
    #     blen[i]=2
      if blen[i] is None:
        if bwidth[i] is not None:
          blen[i]   = int(np.ceil((brange[i][1] + bwidth[i] - brange[i][0])/bwidth[i]))
        else:
          blen[i]   = 10
          bwidth[i] = (brange[i][1] - brange[i][0])/blen[i]
      else:
        bwidth[i] = (brange[i][1] - brange[i][0])/blen[i]
      if wnorm[i] is None: wnorm[i]=bwidth[i]
    #   bins.append(np.linspace(brange[i][0],brange[i][1]+bwidth[i],blen[i]))
      bins.append(np.linspace(brange[i][0],brange[i][1],blen[i]))

    # print(bins)
    # print(brange)
    # print(bwidth)
    # print(blen)
    # Calculate the multi dimensional histo, normalized by wnorm
    h,b=np.histogramdd(axis,weights=w/np.product(wnorm),bins=bins)

    # Return the bins and histo
    return b,h

  def h1(self,axis,bwidth=None,brange=None,wnorm=None,select=None):
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
    hn, h2, h3
    """
    if not brange : brange = [None,None]

    b,h=self.hn([axis],bwidth=[bwidth],brange=[brange],wnorm=[wnorm],select=select)

    # # Verifier b[:-1]
    # h=list(h)
    # h.append(0.0)
    # h=np.array(h)
    return b[0],h

  def h2(self,axis1,axis2,bwidth1=None,bwidth2=None,brange1=None,brange2=None,wnorm=None,select=None):
    """
    Create and return the 2 dimensional histogram of given axis.

    Parameters
    ----------
    axis1, axis2 : str or np.array
      axis to hist
    bwidth1, bwidth2 : float, optional
      bin width. If None, a calculation is done to have 10 bins in the axis
    brange1, brange2 : list of 2 float, optional
      bin maximum and minimum. If a brange element is None, the axis minimum/maximum is taken
    select : dict, optional
      filtering dictionnary

    Returns
    -------
    b1, b2 : np.array
      bins
    h : np.array
      histogram

    Notes
    -----
    the h2 method is just a different way to call the generic method hn

    See Also
    --------
    hn, h1, h3
    """
    if not brange1 : brange1 = [None,None]
    if not brange2 : brange2 = [None,None]

    b,h=self.hn([axis1,axis2],bwidth=[bwidth1,bwidth2],brange=[brange1,brange2],wnorm=wnorm,select=select)

    return b[0],b[1],h

  def h3(self,axis1,axis2,axis3,bwidth1=None,bwidth2=None,bwidth3=None,brange1=None,brange2=None,brange3=None,select=None):
    """
    Create and return the 3 dimensional histogram of given axis.

    Parameters
    ----------
    axis1, axis2, axis3 : str or np.array
      axis to hist
    bwidth1, bwidth2, bwidth3 : float, optional
      bin width. If None, a calculation is done to have 10 bins in the axis
    brange1, brange2, brange3 : list of 2 float, optional
      bin maximum and minimum. If a brange element is None, the axis minimum/maximum is taken
    select : dict, optional
      filtering dictionnary

    Returns
    -------
    b1, b2, b3 : np.array
      bins
    h : np.array
      histogram

    Notes
    -----
    the h3 method is just a different way to call the generic method hn

    See Also
    --------
    hn, h1, h2
    """
    if not brange1 : brange1 = [None,None]
    if not brange2 : brange2 = [None,None]
    if not brange3 : brange3 = [None,None]

    b,h=self.hn([axis1,axis2,axis3],bwidth=[bwidth1,bwidth2,bwidth3],brange=[brange1,brange2,brange3],select=select)

    return b[0],b[1],b[2],h
