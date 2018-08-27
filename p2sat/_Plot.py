#coding:utf8
import numpy as np
import matplotlib.pyplot as plt

class _Plot(object):
  """
  Plots
  """
  def __init__(self,PhaseSpace):
    self._ps=PhaseSpace
    self._h=self._ps.hist
    self.autoclear = False
    self.cmap="viridis"
    plt.ion()

  def get_labels(self,axes,wnorm):
    """
    Returns the labels of given axes.

    Parameters
    ----------
    axes : list of str
      Names of the axes
    wnorm : float or None
      Weight normalization. If None, the last labels element is \"Number\", otherwise it is \"Number/unit1/unit2/...\"

    Returns
    -------
    labels : list of str
      Labels of given axes and label of weight
    """
    names = {
              'x'     : '$x$',
              'y'     : '$y$',
              'z'     : '$z$',
              'px'    : '$p_x$',
              'py'    : '$p_y$',
              'pz'    : '$p_z$',

              'r'     : '$r$',
              'p'     : '$p$',
              'ekin'  : '$E_{kin}$',
              'gamma' : '$\gamma$',

              'theta' : '$\\theta$',
              'phi'   : '$\\phi$'
             }

    units = {
              'x'     : 'um',
              'y'     : 'um',
              'z'     : 'um',
              'px'    : 'MeV/c',
              'py'    : 'MeV/c',
              'pz'    : 'MeV/c',

              'r'     : 'um',
              'p'     : 'MeV/c',
              'ekin'  : 'MeV',
              'gamma' : None,

              'theta' : 'deg',
              'phi'   : 'deg'
              }

    labels=[]
    res=""
    for ax in axes:
      if type(ax) is not str:
        labels.append("")
      else:
        name = self._ps.raw.labels[ax]
        unit = self._ps.raw.units[ax]
        if unit is not None:
          labels.append("{} ({})".format(name,unit))
          if wnorm is None:res += "/{}".format(unit)
        else:
          labels.append(name)

    labels.append("Number{}".format(res))
    return labels

  def figure(self,number=None,clear=True):
    """
    Creates a new figure with given number.

    Parameters
    ----------
    number : int, optional
      Figure number to create
    clear : bool, optional
      Call or not the `clear` method for given number. Default is True

    See Also
    --------
    clear
    """
    plt.figure(number)
    if clear: self.clear(number)

  def clear(self,number=None):
    """
    Clear a plot.

    Parameters
    ----------
    number : int, optional
      Figure number to clear. If None, clear the current figure
    """
    if number is not None:
      plt.figure(number)
    plt.clf()

  def set_title(self,title):
    """

    """
    pass

  def h1(self,axis,where='post',log=False,polar=False,reverse=False,**kargs):
    """
    Plot the 1d histogram of given axis.

    Parameters
    ----------
    axis : str
      Name of the axis to plot
    where : str, optional
      ...
    log : bool, optional
      True to set log scale on y axis
    polar : bool, optional
      True to use a polar plot. axis must be an angle
    reverse : bool, optional
      True to plot axis against number instead of number against axis
    kargs : dict, optional
      Dictionnary to pass to the hist.h1 method

    See Also
    --------
    hist.h1
    """
    if self.autoclear : self.clear()
    if polar:
      a=plt.gca(polar=True)
    else:
      a=plt.gca()

    labels=self.get_labels([axis],kargs.get('wnorm',None))

    b,h=self._h.h1(axis,**kargs)

    if polar: b = np.radians(b)

    if where=='post':
      b=b[:-1]
    elif where=='pre':
      b=b[1:]
    elif where=='mid':
      bsize = b[1:]-b[:-1]
      b = b[:-1] + bsize/2.

    if reverse:
      # Reverse values
      tmp = [b,h]
      h,b = tmp
      # Reverse labels
      tmp=list(labels)
      labels[0]=tmp[1]
      labels[1]=tmp[0]

    a.step(b,h,'.',where=where)

    if polar:
      if log:a.set_rscale('log')
    else:
      a.set_xlim(xmin=min(b),xmax=max(b))
      a.set_xlabel(labels[0])
      a.set_ylabel(labels[1])
      if log:a.set_yscale('log')

    a.grid(True)

    return a

  def f1(self,axis,func_name,log=False,polar=False,reverse=False,**kargs):
    """
    Plot the 1d fit of given axis.

    Parameters
    ----------
    axis : str
      Name of the axis to plot
    func_name : str
      name of the fit function
    log : bool, optional
      True to set log scale on y axis
    polar : bool, optional
      True to use a polar plot. axis must be an angle
    reverse : bool, optional
      True to plot axis against number instead of number against axis
    kargs : dict, optional
      Dictionnary to pass to the hist.h1 method

    See Also
    --------
    hist.f1
    """
    if self.autoclear : self.clear()
    if polar:
      a=plt.gca(polar=True)
    else:
      a=plt.gca()

    labels=self.get_labels([axis],kargs.get('wnorm',None))

    b,h=self._h.f1(axis,func_name,return_fit=True,**kargs)

    if polar: b = np.radians(b)

    if reverse:
      # Reverse values
      tmp = [b,h]
      h,b = tmp
      # Reverse labels
      tmp=list(labels)
      labels[0]=tmp[1]
      labels[1]=tmp[0]

    a.plot(b,h,'-',label="%s fit"%func_name)
    a.legend()
    if polar:
      if log:a.set_rscale('log')
    else:
      a.set_xlim(xmin=min(b),xmax=max(b))
      a.set_xlabel(labels[0])
      a.set_ylabel(labels[1])
      if log:a.set_yscale('log')

    a.grid(True)

    return a

  def h2(self,axis1,axis2,log=False,polar=False,**kargs):
    """
    Plot the 2d histogram of given axes.

    Parameters
    ----------
    axis1,axis2 : str
      Name of the axes to plot
    log : bool, optional
      True to set log scale on y axis
    polar : bool, optional
      True to use a polar plot. axis1 must be an angle
    kargs : dict, optional
      Dictionnary to pass to the hist.h2 method

    See Also
    --------
    hist.h2
    """
    if self.autoclear : self.clear()
    if polar:
      a=plt.gca(polar=True)
    else:
      a=plt.gca()
    labels=self.get_labels([axis1,axis2],kargs.get('wnorm',None))

    b1,b2,h=self._h.h2(axis1,axis2,**kargs)
    g1,g2=np.meshgrid(b1,b2,indexing='ij')

    if log:
      from matplotlib.colors import LogNorm
      norm=LogNorm()
    else:
      norm=None

    if polar:
      a2=a.pcolormesh(np.radians(g1),g2,h,norm=norm,cmap=self.cmap)
    else:
      a2=a.pcolormesh(g1,g2,h,norm=norm,cmap=self.cmap)
      a.set_xlim(xmin=min(b1),xmax=max(b1))
      a.set_xlabel(labels[0])
      a.set_ylabel(labels[1])

    a.grid(True)
    plt.colorbar(a2,label=labels[2])

    return a

  def c2(self,axis1,axis2,log=False,polar=False,gfilter=0.0,**kargs):
    """
    Plot the 2d contour of given axes.

    Parameters
    ----------
    axis1,axis2 : str
      Name of the axes to plot
    log : bool, optional
      True to set log scale on y axis
    polar : bool, optional
      True to use a polar plot. axis1 must be an angle
    gfilter : float, optional
      Filtering scipy.ndimage.filters.gaussian_filter
    kargs : dict, optional
      Dictionnary to pass to the hist.h2 method

    See Also
    --------
    hist.h2
    """
    if self.autoclear : self.clear()
    if polar:
      a=plt.gca(polar=True)
    else:
      a=plt.gca()
    labels=self.get_labels([axis1,axis2],kargs.get('wnorm',None))

    if log:
      from matplotlib.colors import LogNorm
      norm=LogNorm()
    else:
      norm=None

    b1,b2,h=self._h.h2(axis1,axis2,**kargs)

    if polar: b1 = np.radians(b1)

    g1,g2=np.meshgrid(b1,b2,indexing='ij')
    #if gfilter>0.0:
    from scipy.ndimage.filters import gaussian_filter

    a2=a.contour(g1[:-1,:-1],g2[:-1,:-1],gaussian_filter(h,gfilter),
                   norm=norm,colors='k')

    if polar:
      pass
    else:
      a.set_xlim(xmin=min(b1),xmax=max(b1))
      a.set_xlabel(labels[0])
      a.set_ylabel(labels[1])

    plt.clabel(a2, inline=1, fontsize=10 ,fmt='%1.1e')

    return a

  def s2(self,axis1,axis2,log=False,polar=False,select=None):
    """
    Plot the 2d scattering plot of given axes.

    Parameters
    ----------
    axis1,axis2 : str
      Name of the axes to plot
    log : bool, optional
      True to set log scale on y axis
    polar : bool, optional
      True to use a polar plot. axis1 must be an angle
    select : dict, optional
      select dictionnary as in the hist.h2 method
    """
    if self.autoclear : self.clear()
    if polar:
      a=plt.gca(polar=True)
    else:
      a=plt.gca()
    labels=self.get_labels([axis1,axis2],wnorm=1)

    d = self._ps.data

    if type(axis1) is str:axis1 = eval("d.%s"%axis1)
    if type(axis2) is str:axis2 = eval("d.%s"%axis2)
    w   = np.array(d.w)

    if select is not None:
      w = d.select(w,faxis=select.keys(),frange=select.values())
      axis1=d.select(axis1,faxis=select.keys(),frange=select.values())
      axis2=d.select(axis2,faxis=select.keys(),frange=select.values())

    if polar: axis1=np.radians(axis1)

    if log:
      from matplotlib.colors import LogNorm
      norm=LogNorm()
    else:
      norm=None

    a2 = a.scatter(axis1,axis2,c=w,norm=norm,cmap=self.cmap)

    if polar:
      pass
    else:
      a.set_xlabel(labels[0])
      a.set_ylabel(labels[1])

    plt.colorbar(a2,label=labels[2])

    return a

  def h2h1(self,axis1,axis2,log=False,**kargs):
    """
    TODO : doc + kargs + delete labels on h2
    """
    # https://matplotlib.org/examples/pylab_examples/scatter_hist.html
    # https://matplotlib.org/examples/axes_grid/demo_edge_colorbar.html
    kargs1={'X':kargs.get('X',None),
            'erange':kargs.get('erange',[None,None]),
            'bwidth':kargs.get('bwidth1',None),
            'brange':kargs.get('brange1',[None,None])
            }

    kargs2={'X':kargs.get('X',None),
            'erange':kargs.get('erange',[None,None]),
            'bwidth':kargs.get('bwidth2',None),
            'brange':kargs.get('brange2',[None,None])
            }

    tmp = bool(self.autoclear)
    if self.autoclear : self.clear()
    self.autoclear=False

    plt.subplots_adjust(hspace=0.15,wspace=0.15)

    ax1=plt.subplot(221)
    self.h1(axis2,log=False,reverse=True)
    if log:ax1.set_xscale('log')

    ax2=plt.subplot(224)

    self.h1(axis1,log=log)

    ax3=plt.subplot(222,sharex=ax2,sharey=ax1)
    self.h2(axis1,axis2,log=log,**kargs)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)

    self.autoclear=tmp

  def s2h1(self,axis1,axis2,log=False):
    """
    """
    pass

  def h3(self,axis1,axis2,axis3,
        snorm=1.0,hmin=0.0,
        **kargs):
    """
    Plot the 3d histogram of given axes.

    Parameters
    ----------
    axis1,axis2,axis3 : str
      Name of the axes to plot
    kargs : dict, optional
      Dictionnary to pass to the hist.h2 method

    See Also
    --------
    hist.h3
    """
    from mpl_toolkits.mplot3d import Axes3D
    a = plt.subplot(projection='3d')
    b1,b2,b3,h=self._h.h3(axis1,axis2,axis3,**kargs)
    g1,g2,g3=np.meshgrid(b1,b2,b3,indexing='ij')

    tmp=np.array(h)

    for i1,e1 in enumerate(h):
      for i2,e2 in enumerate(e1):
        for i3,e3 in enumerate(e2):
          if e3<hmin:
            tmp[i1][i2][i3]=0.0
          else:
            tmp[i1][i2][i3]=e3

    a.scatter3D(g1,g2,g3,s=snorm*tmp,cmap='hot')
