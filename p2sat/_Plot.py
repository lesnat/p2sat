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
    
  def h1(self,axis, # axis
        label=["",""],log=True, # plot options
        **kargs): # hist options
    """
    
    """
    if self.autoclear : plt.clf()
    b,h=self._h.h1(axis,**kargs)
    plt.step(b,h,'.',label=label[1],where='post') # Verif
    plt.xlim(xmin=min(b),xmax=max(b))
    plt.xlabel(label[0])
    if log:plt.yscale('log')
    plt.legend()

  def h2(self,axis1,axis2,
        label=["","",""],log=False,contour=True,
        **kargs):
    """
    
    """
    if self.autoclear : plt.clf()
    if log:
      from matplotlib.colors import LogNorm
      norm=LogNorm()
    else:
      norm=None
    b1,b2,h=self._h.h2(axis1,axis2,**kargs)
    g1,g2=np.meshgrid(b1,b2,indexing='ij')
    a1=plt.pcolormesh(g1,g2,h,norm=norm) # Voire .T
    
    plt.xlim(xmin=min(b1),xmax=max(b1))
    plt.ylim(ymin=min(b2),ymax=max(b2))
    
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.title(label[2])
    #ca1=plt.colorbar(a1,orientation='horizontal')
    plt.colorbar(a1)
    plt.legend()
    
  def contour(self,axis1,axis2,
              label=["","",""],log=False,gfilter=0.0,
              **kargs):
    """
    
    """
    if self.autoclear : plt.clf()
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
    plt.xlabel(label[0])
    plt.ylabel(label[1])
    plt.title(label[2])
    plt.legend()
    
  
  def h1h2(self,axis1,axis2,
          label=["","",""],log=False,
          **kargs):
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
    if self.autoclear : plt.clf()
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
  
  def scatter(self,axis1,axis2,
              label=["","",""],log=False,
              **kargs):
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()         # no labels
    
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
    
    
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_h2 = [left, bottom, width, height]
    rect_h1x = [left, bottom_h, width, 0.2]
    rect_h1y = [left_h, bottom, 0.2, height]

    axh2 = plt.axes(rect_h2)
    axh1x = plt.axes(rect_h1x)
    axh1y = plt.axes(rect_h1y)

    # no labels
    axh1x.xaxis.set_major_formatter(nullfmt)
    axh1y.yaxis.set_major_formatter(nullfmt)

    # the h2 plot:
    b1,b2,h=self._h.h2(axis1,axis2,**kargs)
    g1,g2=np.meshgrid(b1,b2,indexing='ij')
    axh2.pcolormesh(g1,g2,h) # Voire .T
    ####
    
    """
    # now determine nice limits by hand:
    binwidth = 0.25
    xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
    lim = (int(xymax/binwidth) + 1) * binwidth
    
    axScatter.set_xlim((-lim, lim))
    axScatter.set_ylim((-lim, lim))
    bins = np.arange(-lim, lim + binwidth, binwidth)
    """
    
    ###
    b,h=self._h.h1(axis1,**kargs1)
    axh1x.step(b,h,'.',label=label[1],where='post') # Verif
    b,h=self._h.h1(axis2,**kargs2)
    axh1y.step(h,b,'.',label=label[1],where='post') # Verif
    ####
    """
    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())
    """

  def h3(self,axis1,axis2,axis3,
        snorm=1.0,hmin=0.0,
        **kargs):
    from mpl_toolkits.mplot3d import Axes3D
    ax = plt.subplot(projection='3d')
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
    
    ax.scatter3D(g1,g2,g3,s=snorm*tmp,cmap='hot')

