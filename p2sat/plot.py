#coding:utf8

"""
Plot raw data.

Attributes
----------
autoclear : bool
True to auto clear figure when calling a plot method. Default is False
cmap : str
color map to use in 2d plot. Default is viridris
rcParams : dict
shortcut to matplotlib.rcParams dictionnary (change plot appearence)
"""

from . import hist as _hist

import numpy as _np

import matplotlib.pyplot as _plt
from matplotlib import rcParams

autoclear = False
cmap="viridis"

rcParams = rcParams
rcParams['figure.figsize'] = [9.5, 5.10]
rcParams['mathtext.default']='regular'
rcParams['font.size']=16.0
rcParams['figure.subplot.bottom'] = 0.15

_plt.style.use('bmh')
_plt.ion()

def get_labels(ds, qties, weight, normed):
    r"""
    Returns the labels of given qties.

    Parameters
    ----------
    ds : {PhaseSpace, ScalarField, EventLocation}
        Dataset to use.
    qties : list of str
        Names of the qties
    normed : list of bool or None
        Weight normalization. If None, the last labels element is "Number", otherwise it is "Number/unit1/unit2/..."

    Returns
    -------
    labels : list of str
        Labels of given qty and label of weight
    """
    # Initialization
    labels      = []

    # Get qty labels and units
    qties_labels   = []
    qties_units    = []
    for qty in qties:
        qties_labels.append(ds.metadata.quantity[qty]["label"])
        qty_dim = ds.metadata.quantity[qty]["dimension"]
        qties_units.append(ds.metadata.unit[qty_dim]["label"])

    # Construct axis labels
    axis_labels = []
    for qty_label, qty_unit in zip(qties_labels, qties_units):
        if qty_unit=="":
            axis_labels.append("${label}$".format(label=qty_label))
        else:
            axis_labels.append("${label}$ [${unit}$]".format(label=qty_label,unit=qty_unit))

    # Get weight label and unit
    weight_label = ds.metadata.quantity[weight]["label"]
    weight_dim   = ds.metadata.quantity[weight]["dimension"]
    weight_unit  = ds.metadata.unit[weight_dim]["label"]

    # Construct "title" label
    title_unit  = " ~".join(qties_units)
    if normed == False:
        title_label = "${label}$".format(label=weight_label)
    else:
        numerator_label = "d" + weight_label
        denominator_label = "d" + " ~d".join(qties_labels)
        if title_unit.strip() == "":
            title_label = "$%s/%s$"%(numerator_label, denominator_label)
        else:
            title_label = "$%s/%s$ [$(%s)^{-1}$]"%(numerator_label, denominator_label, title_unit)

    return [*axis_labels, title_label]

def figure(number=None,clear=False):
    r"""
    Creates a new figure with given number.

    Parameters
    ----------
    number : int, optional
        Figure number to create
    clear : bool, optional
        Call or not the `clear` method for given number. Default is True

    See Also
    --------
    plot.clear
    """
    _plt.figure(number)
    if clear: clear_figure(number)

def clear_figure(number=None):
    r"""
    Clear a plot.

    Parameters
    ----------
    number : int, optional
        Figure number to clear. If None, clear the current figure
    """
    if number is not None:
        _plt.figure(number)
    _plt.clf()

def set_title(title,number=None):
    r"""
    Set the title of the figure.
    """
    if number is not None:
        figure(number,clear=False)
    _plt.title(title)

def hist1d(ds, qty, where='post', legend="", log=False, polar=False, reverse=False, clear=False,**kargs):
    r"""
    Plot the 1d histogram of given qty.

    Parameters
    ----------
    ds : {PhaseSpace, ScalarField, EventLocation}
        Dataset to use.
    qty : str
        Name of the qty to plot
    where : str, optional
        ...
    log : bool, optional
        True to set log scale on y qty
    polar : bool, optional
        True to use a polar plot. qty must be an angle
    reverse : bool, optional
        True to plot qty against number instead of number against qty
    clear: bool, optional
        Clear or not the figure before plotting
    kargs : dict, optional
        Dictionnary to pass to the hist.hist1d method

    See Also
    --------
    hist.hist1d
    """
    if autoclear or clear: clear_figure()
    if polar:
        a=_plt.gca(polar=True)
    else:
        a=_plt.gca()

    labels=get_labels(ds,[qty],kargs.get('weight','w'),kargs.get('normed',True))

    b,h=_hist.hist1d(ds, qty, **kargs)

    if polar: b = _np.radians(b) # Polar plot usefull ?????????

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

    if legend != "":
        a.step(b,h,'.-',where=where, label=legend)
        _plt.legend()
    else:
        a.step(b,h,'.-',where=where)

    if polar:
        if log:a.set_rscale('log')
    else:
        a.set_xlim(xmin=min(b),xmax=max(b))
        a.set_xlabel(labels[0])
        a.set_ylabel(labels[1])
        if log:a.set_yscale('log')

    a.grid(True)

    _plt.show()

    return a

# def fit1d(ds, qty, func_name, log=False, polar=False, reverse=False, clear=False,**kargs):
#     r"""
#     Plot the 1d fit of given qty.
#
#     Parameters
#     ----------
#     ds : {PhaseSpace, ScalarField, EventLocation}
#         Dataset to use.
#     qty : str
#         Name of the qty to plot
#     func_name : str
#         name of the fit function
#     log : bool, optional
#         True to set log scale on y qty
#     polar : bool, optional
#         True to use a polar plot. qty must be an angle
#     reverse : bool, optional
#         True to plot qty against number instead of number against qty
#     clear: bool, optional
#         Clear or not the figure before plotting
#     kargs : dict, optional
#         Dictionnary to pass to the hist.hist1d method
#
#     See Also
#     --------
#     hist.fit1d
#     """
#     if autoclear or clear: clear_figure()
#     if polar:
#         a=_plt.gca(polar=True)
#     else:
#         a=_plt.gca()
#
#     labels=get_labels(ds, [qty],kargs.get('weight','w'),kargs.get('normed',True))
#
#     b,h=_hist.fit1d(qty,func_name,return_fit=True,**kargs)
#
#     if polar: b = _np.radians(b)
#
#     if reverse:
#         # Reverse values
#         tmp = [b,h]
#         h,b = tmp
#         # Reverse labels
#         tmp=list(labels)
#         labels[0]=tmp[1]
#         labels[1]=tmp[0]
#
#     a.plot(b,h,'-',label="%s fit"%func_name)
#     a.legend()
#     if polar:
#         if log:a.set_rscale('log')
#     else:
#         #a.set_xlim(xmin=min(b),xmax=max(b))
#         a.set_xlabel(labels[0])
#         a.set_ylabel(labels[1])
#         if log:a.set_yscale('log')
#
#     a.grid(True)
#
#     _plt.show()
#
#     return a

# def animate1d(ds, qty, t=None, pause=.5, where='post', log=False, polar=False, reverse=False,**kargs):
#     r"""
#     Plot the animation of 1d histogram of given qty.
#
#     Parameters
#     ----------
#     ds : {PhaseSpace, ScalarField, EventLocation}
#         Dataset to use.
#     qty : str
#         Name of the qty to plot
#     t : dict
#         parameters to construct plot times
#     pause : float
#         waiting time between each plot (in seconds)
#     where : str, optional
#         ...
#     log : bool, optional
#         True to set log scale on y qty
#     polar : bool, optional
#         True to use a polar plot. qty must be an angle
#     reverse : bool, optional
#         True to plot qty against number instead of number against qty
#     kargs : dict, optional
#         Dictionnary to pass to the hist.hist1d method
#
#     Notes
#     -----
#     t dictionnary must contains following keys :
#
#     - "min", to define the first plot time
#     - "max", to define the last plot time
#     - "Nbins", to define the number of timesteps to plot
#
#     Examples
#     --------
#     >>> eps = ExamplePhaseSpace()
#     >>> eps.plot.autoclear = False
#     >>> #eps.plot.animate1d('r',t={'min':0,'max':1000,'Nbins':10},bwidth=1.)
#
#     See Also
#     --------
#     plot.hist1d
#     """
#     copy = ds._ps.copy()
#     a = _plt.gca()
#
#     T = _np.linspace(t["min"],t["max"],t["Nbins"])
#
#     for tt in T:
#         copy.data.propagate(t=tt,verbose=True)
#         if autoclear:a.clear_figure()
#         a=hist1d(copy, qty,where=where,log=log,polar=polar,reverse=reverse,**kargs)
#         a.set_label('t = %.3E fs'%tt)
#         _plt.pause(pause)
#     _plt.legend()
#
#     return a

def hist2d(ds, qty1, qty2, log=False, polar=False, clear=False, **kargs):
    r"""
    Plot the 2d histogram of given qty.

    Parameters
    ----------
    ds : {PhaseSpace, ScalarField, EventLocation}
        Dataset to use.
    qty1,qty2 : str
        Name of the qty to plot
    log : bool, optional
        True to set log scale on y qty
    polar : bool, optional
        True to use a polar plot. qty1 must be an angle
    clear : bool, optional
        Clear or not the figure before plotting
    kargs : dict, optional
        Dictionnary to pass to the hist.hist2d method

    See Also
    --------
    hist.hist2d
    """
    if autoclear or clear: clear_figure()
    if polar:
        a=_plt.gca(polar=True)
    else:
        a=_plt.gca()
    labels=get_labels(ds, [qty1,qty2],kargs.get('weight','w'),kargs.get('normed',True))

    b1,b2,h=_hist.hist2d(ds,qty1,qty2,**kargs)
    g1,g2=_np.meshgrid(b1,b2,indexing='ij')

    if log:
        from matplotlib.colors import LogNorm
        norm=LogNorm(vmax=h.max())
    else:
        from matplotlib.colors import Normalize
        norm=Normalize(vmax=_np.max(h))
        norm=None

    if polar:
        a2=a.pcolormesh(_np.radians(g1),g2,h,norm=norm,cmap=cmap)
        _plt.colorbar(a2,label=labels[2])
    else:
        a2=a.pcolormesh(g1,g2,h,norm=norm,cmap=cmap)
        a.set_xlim(xmin=min(b1),xmax=max(b1))
        a.set_ylim(ymin=min(b2),ymax=max(b2))
        a.set_xlabel(labels[0])
        a.set_ylabel(labels[1])
        _plt.colorbar(a2,label=labels[2])

    a.grid(True)

    _plt.show()

    return a

# def trisurf2d(ds,qty1,qty2,log=False,polar=False,clear=False,spherical=False,**kargs):
#     r"""
#     Plot the 2d triangulated surface of given qty.
#
#     Parameters
#     ----------
#     ds : {PhaseSpace, ScalarField, EventLocation}
#         Dataset to use.
#     qty1,qty2 : str
#         Name of the qty to plot
#     log : bool, optional
#         True to set log scale on y qty
#     polar : bool, optional
#         True to use a polar plot. qty1 must be an angle
#     clear : bool, optional
#         Clear or not the figure before plotting
#     spherical : bool, optional
#         Plot the projection on the unit sphere. Must be done with (qty1, qty2) = (theta,phi)
#     kargs : dict, optional
#         Dictionnary to pass to the hist.hist2d method
#
#     See Also
#     --------
#     hist.hist2d
#     """
#     if autoclear or clear: clear_figure()
#
#     from mpl_toolkits.mplot3d import Axes3D
#     a=_plt.gca(projection='3d')
#
#     labels=get_labels(ds, [qty1,qty2],kargs.get('weight','w'),kargs.get('normed',True))
#
#     b1,b2,h=_hist.hist2d(ds,qty1,qty2,**kargs)
#     g1,g2=_np.meshgrid(b1,b2,indexing='ij')
#
#     if log:
#         from matplotlib.colors import LogNorm
#         norm=LogNorm(vmax=h.max())
#     else:
#         from matplotlib.colors import Normalize
#         norm=Normalize(vmax=_np.max(h))
#         norm=None
#
#     # https://stackoverflow.com/questions/6539944/color-matplotlib-plot-surface-command-with-surface-gradient/6543777#6543777
#     from matplotlib import cm
#     # r = ds.read
#     #
#     # btheta = 1./_np.degrees(_np.sin(list(reversed(_np.linspace(1e-5,_np.pi/2.-1e-5,100)))))
#     # bphi = _np.linspace(1e-5,2*180-1e-5,100)
#     #
#     # h, theta, phi = _np.histogram2d(r.theta, r.phi, weights=r.w, bins=(btheta,bphi))
#     # g1,g2=_np.meshgrid(theta,phi,indexing='ij')
#
#     g1 = _np.radians(g1)
#     g2 = _np.radians(g2)
#     h = h/_np.sin(g1[1:,1:])
#     # filtr = _np.nonzero(h)
#     # h = h[filtr]
#     # Delete extreme values
#     # g1 = _np.radians(g1[1:-1,1:-1])
#     # g2 = _np.radians(g2[1:-1,1:-1])
#     # h = h[1:,1:]/_np.sin(g1)
#     # # We choose to delete theta=180
#     # h = h[:-1,:-1]
#
#     x = _np.sin(g1) * _np.cos(g2)
#     y = _np.sin(g1) * _np.sin(g2)
#     z = _np.cos(g1)
#
#     # a2 = a.plot_trisurf(x[filtr], y[filtr], z[filtr], facecolors=cm.viridis(h/h.max()), color='c', alpha=0.6, norm=norm, linewidth=0)
#     a2 = a.contourf(x, y, z, facecolors=cm.viridis(h/h.max()), color='c', alpha=0.6, norm=norm, linewidth=0)
#     # a2 = a.contour3D(x, y, z, facecolors=cm.viridis(h/h.max()), color='c', alpha=0.6, norm=norm, linewidth=0)
#
#     a.set_xlim([-1,1])
#     a.set_ylim([-1,1])
#     a.set_zlim([-1,1])
#     m = cm.ScalarMappable(cmap=_plt.cm.viridis, norm=norm)
#     m.set_array([])
#     _plt.colorbar(m)
#     a.set_xlabel(r"$\sin \theta \cos \phi$")
#     a.set_ylabel(r"$\sin \theta \sin \phi$")
#     a.set_zlabel(r"$\cos \theta$")
#
#     a.grid(True)
#
#     _plt.show()
#
#     return a

# def animate2d(ds,qty1,qty2,t=None,pause=.5,log=False,polar=False,**kargs):
#     r"""
#     Plot the animation of 2d histogram of given qty.
#
#     Parameters
#     ----------
#     ds : {PhaseSpace, ScalarField, EventLocation}
#         Dataset to use.
#     qty1,qty2 : str
#         Name of the qty to plot
#     t : dict
#         parameters to construct plot times
#     pause : float
#         waiting time between each plot (in seconds)
#     log : bool, optional
#         True to set log scale on y qty
#     polar : bool, optional
#         True to use a polar plot. qty1 must be an angle
#     kargs : dict, optional
#         Dictionnary to pass to the hist.hist2d method
#
#     Notes
#     -----
#     t dictionnary must contains following keys :
#
#     - "min", to define the first plot time
#     - "max", to define the last plot time
#     - "Nbins", to define the number of timesteps to plot
#
#     Examples
#     --------
#
#     >>> eps.plot.a2('x','y',t={'min':0,'max':1000,'Nbins':10},log=True) #doctest: +SKIP
#
#     See Also
#     --------
#     plot.hist2d
#     """
#     copy = ds._ps.copy()
#     a = _plt.gca()
#
#     T = _np.linspace(t["min"],t["max"],t["Nbins"])
#
#     for tt in T:
#         copy.data.propagate(t=tt,verbose=True)
#         _plt.clf()
#         a=copy.plot.hist2d(qty1,qty2,log=log,polar=polar,**kargs)
#         a.set_title('t = %.3E fs'%tt)
#         _plt.pause(pause)
#
#     return a

# def contour2d(ds,qty1,qty2,log=False,polar=False,gfilter=0.0,clear=False,**kargs):
#     r"""
#     Plot the 2d contour of given qty.
#
#     Parameters
#     ----------
#     ds : {PhaseSpace, ScalarField, EventLocation}
#         Dataset to use.
#     qty1,qty2 : str
#         Name of the qty to plot
#     log : bool, optional
#         True to set log scale on y qty
#     polar : bool, optional
#         True to use a polar plot. qty1 must be an angle
#     gfilter : float, optional
#         Filtering scipy.ndimage.filters.gaussian_filter
#     clear: bool, optional
#         Clear or not the figure before plotting
#     kargs : dict, optional
#         Dictionnary to pass to the hist.hist2d method
#
#     See Also
#     --------
#     hist.hist2d
#     """
#     if autoclear or clear: clear_figure()
#     if polar:
#         a=_plt.gca(polar=True)
#     else:
#         a=_plt.gca()
#     labels=get_labels(ds, [qty1,qty2],kargs.get('weight','w'),kargs.get('normed',True))
#
#     if log:
#         from matplotlib.colors import LogNorm
#         norm=LogNorm()
#     else:
#         norm=None
#
#     b1,b2,h=_hist.hist2d(ds,qty1,qty2,**kargs)
#
#     if polar: b1 = _np.radians(b1)
#
#     g1,g2=_np.meshgrid(b1,b2,indexing='ij')
#     #if gfilter>0.0:
#     from scipy.ndimage.filters import gaussian_filter
#
#     a2=a.contour(g1[:-1,:-1],g2[:-1,:-1],gaussian_filter(h,gfilter),
#                 norm=norm,colors='k')
#
#     if polar:
#         pass
#     else:
#         a.set_xlim(xmin=min(b1),xmax=max(b1))
#         a.set_xlabel(labels[0])
#         a.set_ylabel(labels[1])
#
#     _plt.clabel(a2, inline=1, fontsize=10 ,fmt='%1.1e')
#
#     _plt.show()
#
#     return a

# def scatter2d(ds,qty1,qty2,weight='w',snorm=1.,log=False,polar=False,clear=False,select=None):
#     r"""
#     Plot the 2d scattering plot of given qty.
#
#     Parameters
#     ----------
#     ds : {PhaseSpace, ScalarField, EventLocation}
#         Dataset to use.
#     qty1,qty2 : str
#         name of the qty to plot
#     weight : str, optional
#         weight to plot. Default is w
#     snorm : float, optional
#         dots size normalization. Default is 1
#     log : bool, optional
#         True to set log scale on y qty
#     polar : bool, optional
#         True to use a polar plot. qty1 must be an angle
#     clear: bool, optional
#         Clear or not the figure before plotting
#     select : dict, optional
#         select dictionnary as in the hist.hist2d method
#     """
#     if autoclear or clear: clear_figure()
#     if polar:
#         a=_plt.gca(polar=True)
#     else:
#         a=_plt.gca()
#     labels=get_labels(ds, [qty1,qty2],weight,normed=False)
#
#     d = ds._ps.data
#
#     qty1   = d.quantity(qty1,select=select)
#     qty2   = d.quantity(qty2,select=select)
#     w       = d.quantity(weight,select=select)
#
#     if polar: qty1=_np.radians(qty1)
#
#     if log:
#         from matplotlib.colors import LogNorm
#         norm=LogNorm()
#     else:
#         norm=None
#
#     s = snorm * w/max(w)
#
#     a2 = a.scatter(qty1,qty2,c=w,s=s,norm=norm,cmap=ds.cmap)
#
#     if polar:
#         pass
#     else:
#         a.set_xlabel(labels[0])
#         a.set_ylabel(labels[1])
#
#     _plt.colorbar(a2,label=labels[2])
#
#     _plt.show()
#
#     return a

# def angles(ds, btheta=5., bphi=5., log=False):
#     a=_plt.gca()
#     labels=get_labels(ds, ["theta","phi"],"w",normed=True)
#
#     bintheta = _np.linspace(0,180,int(180./btheta))
#     binphi = _np.linspace(0,360,int(360./bphi))
#
#     r = ds.read
#     h, bt, bp = _np.histogram2d(r.theta, r.phi, weights=r.w/_np.sin(r.theta), bins=(bintheta,binphi))
#
#     if log:
#         from matplotlib.colors import LogNorm
#         norm = LogNorm()
#     else:
#         norm = None
#
#     gtheta, gphi = _np.meshgrid(bt, bp, indexing='ij')
#     a2 = _plt.pcolormesh(gtheta, gphi, h, norm=norm)
#     a.set_xlabel(labels[0])
#     a.set_ylabel(labels[1])
#     _plt.colorbar(a2,label=r"$d N / d\Omega$ [a.u.]")
#     _plt.show()

# def h2h1(ds,qty1,qty2,log=False,clear=False,**kargs):
#     """
#     TODO
#     """
#     # # https://matplotlib.org/examples/pylab_examples/scatter_hist.html
#     # # https://matplotlib.org/examples/axes_grid/demo_edge_colorbar.html
#     # kargs1={'X':kargs.get('X',None),
#     #         'erange':kargs.get('erange',[None,None]),
#     #         'bwidth':kargs.get('bwidth1',None),
#     #         'brange':kargs.get('brange1',[None,None])
#     #         }
#     #
#     # kargs2={'X':kargs.get('X',None),
#     #         'erange':kargs.get('erange',[None,None]),
#     #         'bwidth':kargs.get('bwidth2',None),
#     #         'brange':kargs.get('brange2',[None,None])
#     #         }
#     #
#     # tmp = bool(autoclear)
#     # if autoclear : clear_figure()
#     # autoclear=False
#     #
#     # _plt.subplots_adjust(hspace=0.15,wspace=0.15)
#     #
#     # ax1=_plt.subplot(221)
#     # ds.hist1d(qty2,log=False,reverse=True)
#     # if log:ax1.set_xscale('log')
#     #
#     # ax2=_plt.subplot(224)
#     #
#     # ds.hist1d(qty1,log=log)
#     #
#     # ax3=_plt.subplot(222,sharex=ax2,sharey=ax1)
#     # ds.hist2d(qty1,qty2,log=log,**kargs)
#     # _plt.setp(ax3.get_yticklabels(), visible=False)
#     # _plt.setp(ax3.get_xticklabels(), visible=False)
#     #
#     # autoclear=tmp
#     #
#     # _plt.show()
#     pass

# def s2h1(ds,qty1,qty2,log=False,clear=False):
#     """
#     TODO
#     """
#     pass

# def hist3d(ds,qty1,qty2,qty3,s=5,wmin=0,log=False,clear=False,**kargs):
#     r"""
#     Plot the 3d histogram of given qty.
#
#     Parameters
#     ----------
#     ds : {PhaseSpace, ScalarField, EventLocation}
#         Dataset to use.
#     qty1,qty2,qty3 : str
#         Name of the qty to plot
#     s : float, optional
#         square sizes. Default is 5
#     wmin : float, optional
#         minimum weight to plot. Default is 0
#     log : bool, optional
#         log color scale. Default is False
#     clear: bool, optional
#         Clear or not the figure before plotting
#     kargs : dict, optional
#         Dictionnary to pass to the hist.hist3d method
#
#     See Also
#     --------
#     hist.hist3d
#     """
#     if autoclear or clear: clear_figure()
#     from mpl_toolkits.mplot3d import Axes3D
#     # https://codereview.stackexchange.com/questions/62180/colorbar-for-matplotlib-3d-patch-plot
#     from matplotlib import cm
#
#     a = _plt.subplot(projection='3d')
#
#     labels=get_labels(ds, [qty1,qty2,qty3],kargs.get('weight','w'),kargs.get('normed',True))
#
#     b1,b2,b3,h=_hist.hist3d(qty1,qty2,qty3,**kargs)
#
#     xs,ys,zs,w = [],[],[],[]
#
#     for i1,e1 in enumerate(h):
#         for i2,e2 in enumerate(e1):
#             for i3,e3 in enumerate(e2):
#                 if e3>wmin:
#                     xs.append(b1[i1])
#                     ys.append(b2[i2])
#                     zs.append(b3[i3])
#                     w.append(h[i1][i2][i3])
#
#     if log:
#         c = _np.log10(w)
#     else:
#         c = w
#
#     a.scatter(xs=xs,ys=ys,zs=zs,c=c,s=s,marker='s')
#
#     a.set_xlim(min(xs),max(xs))
#     a.set_ylim(min(ys),max(ys))
#     a.set_zlim(min(zs),max(zs))
#     a.set_xlabel(labels[0],labelpad=15)
#     a.set_ylabel(labels[1],labelpad=15)
#     a.set_zlabel(labels[2],labelpad=15)
#
#     m = cm.ScalarMappable(cmap=ds.cmap)
#     m.set_array([min(w),max(w)])
#
#     if log:
#         from matplotlib.colors import LogNorm
#         norm=LogNorm()
#     else:
#         norm=None
#
#     m.set_norm(norm)
#
#     f=_plt.gcf()
#
#     f.colorbar(m,label=labels[3])
#
#     _plt.show()
#     return a


# def scatter3d(ds,qty1,qty2,qty3,weight="w",snorm=1.0,log=False,clear=False,select=None):
#     r"""
#     Plot the 3d scatter plot of given qty.
#
#     Parameters
#     ----------
#     ds : {PhaseSpace, ScalarField, EventLocation}
#         Dataset to use.
#     qty1,qty2,qty3 : str
#         Name of the qty to plot
#     weight : str, optional
#         weight to plot. Default is w
#     snorm : float, optional
#         dots size normalization. Default is 1
#     log : bool, optional
#         log color scale. Default is False
#     clear: bool, optional
#         Clear or not the figure before plotting
#     select : dict, optional
#         filtering dictionnary
#     """
#     if autoclear or clear: clear_figure()
#     from mpl_toolkits.mplot3d import Axes3D
#     from matplotlib import cm
#
#     a = _plt.subplot(projection='3d')
#
#     labels=get_labels(ds, [qty1,qty2,qty3],weight,normed=False)
#
#     b1  = ds.read.quantity(qty1,select=select)
#     b2  = ds.read.quantity(qty2,select=select)
#     b3  = ds.read.quantity(qty3,select=select)
#     w   = ds.read.quantity(weight,select=select)
#
#     if log:
#         c = _np.log10(w)
#     else:
#         c = w
#     s = snorm * w/max(w)
#
#     a.scatter(xs=b1,ys=b2,zs=b3,c=c,s=s,marker='o')
#
#     a.set_xlim(min(b1),max(b1))
#     a.set_ylim(min(b2),max(b2))
#     a.set_zlim(min(b3),max(b3))
#     a.set_xlabel(labels[0],labelpad=15)
#     a.set_ylabel(labels[1],labelpad=15)
#     a.set_zlabel(labels[2],labelpad=15)
#
#     m = cm.ScalarMappable(cmap=ds.cmap)
#     m.set_array([min(w),max(w)])
#
#     if log:
#         from matplotlib.colors import LogNorm
#         norm=LogNorm()
#     else:
#         norm=None
#
#     m.set_norm(norm)
#
#     f=_plt.gcf()
#
#     f.colorbar(m,label=labels[3])
#
#     _plt.show()
#     return a
