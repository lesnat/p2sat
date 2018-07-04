#coding:utf8
import numpy as np
import matplotlib.pyplot as plt

class _Plot(object):
  """
  Plots
  """
  def __init__(self,PhaseSpace):
    self._ps=PhaseSpace
    self._r=self._ps.raw
    self._h=self._ps.hist
    self.autoclear = True

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
              'x'     : 'µm',
              'y'     : 'µm',
              'z'     : 'µm',
              'px'    : 'MeV/c',
              'py'    : 'MeV/c',
              'pz'    : 'MeV/c',

              'r'     : 'µm',
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
        name = names[ax]
        unit = units[ax]
        if unit is not None:
          labels.append("{} ({})".format(name,unit))
          if wnorm is None:res += "/{}".format(unit)
        else:
          labels.append(name)

    labels.append("Number{}".format(res))
    return labels

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

  def h1(self,axis,log=True,polar=False,**kargs):
    """
    Plot the 1d histogram of given axis.

    Parameters
    ----------
    axis : str
      Name of the axis to plot
    log : bool, optional
      True to set log scale on y axis
    kargs : dict, optional
      Dictionnary to pass to the hist.h1 method

    See Also
    --------
    hist.h1
    """
    if self.autoclear : self.clear()
    if polar:
      a=plt.axes(polar=True)
    else:
      a=plt.axes()

    labels=self.get_labels([axis],kargs.get('wnorm',None))

    b,h=self._h.h1(axis,**kargs)

    if polar:
      a.step(b[:-1],h,'.-') # Verif
    #   a.step(h,b[:-1],'.-') # Verif
    #   a.set_thetalim(thetamin=min(b),thetamax=max(b))
    #   a.set_rgrids()
    #   a.set_rlabel(labels[0])
    #   if log:a.set_rscale('log')
    else:
      a.step(b[:-1],h,'.',where='post') # Verif
      a.set_xlim(xmin=min(b),xmax=max(b))
      a.set_xlabel(labels[0])
      a.set_ylabel(labels[1])
      if log:a.set_yscale('log')
    #   a.set_grid()

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
    kargs : dict, optional
      Dictionnary to pass to the hist.h2 method

    See Also
    --------
    hist.h2
    """
    if self.autoclear : self.clear()
    if polar:plt.axes(polar=True)
    labels=self.get_labels([axis1,axis2],kargs.get('wnorm',None))

    b1,b2,h=self._h.h2(axis1,axis2,**kargs)
    g1,g2=np.meshgrid(b1,b2,indexing='ij')

    if log:
      from matplotlib.colors import LogNorm
      norm=LogNorm()
    else:
      norm=None

    a1=plt.pcolormesh(g1,g2,h,norm=norm)

    if not polar:
      plt.xlim(xmin=min(b1),xmax=max(b1))
      plt.ylim(ymin=min(b2),ymax=max(b2))
      plt.grid()
      plt.xlabel(labels[0])
      plt.ylabel(labels[1])
      plt.colorbar(a1,label=labels[2])
    else:
      pass
    #   plt.rgrids()


  def c2(self,axis1,axis2,log=False,polar=False,gfilter=0.0,**kargs):
    """
    Plot the 2d contour of given axes.

    Parameters
    ----------
    axis1,axis2 : str
      Name of the axes to plot
    log : bool, optional
      True to set log scale on y axis
    gfilter : float, optional
      Filtering scipy.ndimage.filters.gaussian_filter
    kargs : dict, optional
      Dictionnary to pass to the hist.h2 method

    See Also
    --------
    hist.h2
    """
    if self.autoclear : self.clear()
    if polar:plt.subplot(111,projection='polar')
    labels=self.get_labels([axis1,axis2],kargs.get('wnorm',None))

    if log:
      from matplotlib.colors import LogNorm
      norm=LogNorm()
    else:
      norm=None
    b1,b2,h=self._h.h2(axis1,axis2,**kargs)
    g1,g2=np.meshgrid(b1,b2,indexing='ij')
    #if gfilter>0.0:
    from scipy.ndimage.filters import gaussian_filter
    a2=plt.contour(g1[:-1,:-1],g2[:-1,:-1],gaussian_filter(h,gfilter),
                   norm=norm,colors='k')
    plt.clabel(a2, inline=1, fontsize=10 ,fmt='%1.1e')

    plt.xlim(xmin=min(b1),xmax=max(b1))
    plt.ylim(ymin=min(b2),ymax=max(b2))
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    #plt.title(labels[2])
    plt.legend()

  def s2(self,axis1,axis2,log=False,polar=False,select=None):
    """
    Plot the 2d scattering plot of given axes.

    Parameters
    ----------
    axis1,axis2 : str
      Name of the axes to plot
    """
    if self.autoclear : self.clear()
    if polar:plt.subplot(111,projection='polar')
    r = self._ps.raw
    labels=self.get_labels([axis1,axis2],wnorm=1)
    if type(axis1) is str:axis1 = eval("r.%s"%axis1)
    if type(axis2) is str:axis2 = eval("r.%s"%axis2)
    w   = np.array(r.w)

    if select is not None:
      w = r.select(w,faxis=select.keys(),frange=select.values())
      axis1=r.select(axis1,faxis=select.keys(),frange=select.values())
      axis2=r.select(axis2,faxis=select.keys(),frange=select.values())
    if log:
      from matplotlib.colors import LogNorm
      norm=LogNorm()
    else:
      norm=None
    plt.scatter(axis1,axis2,c=w,norm=norm)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.colorbar(label=labels[2])

  def h2h1(self,axis1,axis2,log=False,**kargs):
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
    b,h=self._h.h1(axis1)
    plt.step(h,b,label=label[2],where='post') # Verif
    plt.ylabel(label[0])
    if log:plt.xscale('log')
    plt.legend()

    ax2=plt.subplot(224)
    self.h1(axis2,label=[label[1],label[2]],log=log)

    ax3=plt.subplot(222,sharex=ax2,sharey=ax1)
    self.h2(axis1,axis2,label=label,log=log,**kargs)
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
