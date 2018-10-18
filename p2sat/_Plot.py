#coding:utf8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

class _Plot(object):
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
    def __init__(self,PhaseSpace):
        self._ps=PhaseSpace
        self._h=self._ps.hist
        self.autoclear = False
        self.cmap="viridis"
        self.rcParams = rcParams
        plt.style.use('bmh')
        self.rcParams['figure.figsize'] = [9.5, 5.10]
        self.rcParams['mathtext.default']='regular'
        self.rcParams['font.size']=16.0
        self.rcParams['figure.subplot.bottom'] = 0.15
        plt.ion()

    def get_labels(self,axes,weight,normed):
        """
        Returns the labels of given axes.

        Parameters
        ----------
        axes : list of str
            Names of the axes
        normed : list of bool or None
            Weight normalization. If None, the last labels element is \"Number\", otherwise it is \"Number/unit1/unit2/...\"

        Returns
        -------
        labels : list of str
            Labels of given axes and label of weight
        """
        # Initialization
        r = self._ps.data.raw
        labels  = []
        ax_label=[]
        ax_unit =[]
        w_label =[]
        w_unit  =[]

        # Get axes labels and units
        for ax in axes:
            if type(ax) is not str:
                ax_label.append("")
                ax_unit.append("")
            else:
                ax_label.append(r.labels[ax])
                ax_unit.append(r.units[ax])

        # Get normalization bool list
        if type(normed) is bool:
            if normed:
                normed=[True]*len(axes)
            else:
                normed=[False]*len(axes)

        # Construct denominator of weight label and unit
        den_label = ""
        den_unit = ""
        for i,ax in enumerate(axes):
            if normed[i]:
                den_label += "d"+ax_label[i]+" \ "
                den_unit += ax_unit[i]+" \ "
        # Delete the last spaces if needed
        if den_label !="":
            den_unit = den_unit[:-3]
            # den_label = den_label[:-3]

        # Construct weight label
        if den_label == "":
            w_label = "$%s$"%r.labels[weight]
        else:
            w_label = "$d%s/%s$"%(r.labels[weight],den_label)

        # Construct weight unit
        if den_unit == "" and r.units[weight] is None:
            w_unit  = ""

        if den_unit == "" and r.units[weight] is not None:
            w_unit = "$%s$"%r.units[weight]

        if den_unit != "" and r.units[weight] is None:
            w_unit  = "$(%s)^{-1}$"%(den_unit)

        if den_unit != "" and r.units[weight] is not None:
            w_unit  = "$(%s)/(%s)$"%(r.units[weight],den_unit)

        # Construct and return labels list
        for i,_ in enumerate(axes):
            labels.append("$%s$ [$%s$]"%(ax_label[i],ax_unit[i]))

        if w_unit == "":
            labels.append(w_label)
        else:
            labels.append("%s [%s]"%(w_label,w_unit))

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
        plot.clear
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

    def set_title(self,title,number=None):
        """
        Set the title of the figure.
        """
        if number is not None:
            self.figure(number,clear=False)
        plt.title(title)

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

        labels=self.get_labels([axis],kargs.get('weight','w'),kargs.get('normed',True))

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

        plt.show()

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

        labels=self.get_labels([axis],kargs.get('weight','w'),kargs.get('normed',True))

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

        plt.show()

        return a

    def a1(self,axis,t=None,pause=.5,where='post',log=False,polar=False,reverse=False,**kargs):
        """
        Plot the animation of 1d histogram of given axis.

        Parameters
        ----------
        axis : str
            Name of the axis to plot
        t : dict
            parameters to construct plot times
        pause : float
            waiting time between each plot (in seconds)
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

        Notes
        -----
        t dictionnary must contains following keys :

        - "min", to define the first plot time
        - "max", to define the last plot time
        - "Nbins", to define the number of timesteps to plot

        Examples
        --------

        >>> eps.plot.autoclear = False
        >>> eps.plot.a1('r',t={'min':0,'max':1000,'Nbins':10},bwidth=1.)

        See Also
        --------
        plot.h1
        """
        copy = self._ps.copy()
        a = plt.gca()

        T = np.linspace(t["min"],t["max"],t["Nbins"])

        for tt in T:
            copy.data.propagate(t=tt,verbose=True)
            if self.autoclear:a.clear()
            a=copy.plot.h1(axis,where=where,log=log,polar=polar,reverse=reverse,**kargs)
            a.set_label('t = %.3E fs'%tt)
            plt.pause(pause)
        plt.legend()

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
        labels=self.get_labels([axis1,axis2],kargs.get('weight','w'),kargs.get('normed',True))

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
            a.set_ylim(ymin=min(b2),ymax=max(b2))
            a.set_xlabel(labels[0])
            a.set_ylabel(labels[1])

        a.grid(True)
        plt.colorbar(a2,label=labels[2])

        plt.show()

        return a

    def a2(self,axis1,axis2,t=None,pause=.5,log=False,polar=False,**kargs):
        """
        Plot the animation of 2d histogram of given axes.

        Parameters
        ----------
        axis1,axis2 : str
            Name of the axes to plot
        t : dict
            parameters to construct plot times
        pause : float
            waiting time between each plot (in seconds)
        log : bool, optional
            True to set log scale on y axis
        polar : bool, optional
            True to use a polar plot. axis1 must be an angle
        kargs : dict, optional
            Dictionnary to pass to the hist.h2 method

        Notes
        -----
        t dictionnary must contains following keys :

        - "min", to define the first plot time
        - "max", to define the last plot time
        - "Nbins", to define the number of timesteps to plot

        Examples
        --------

        >>> eps.plot.a2('x','y',t={'min':0,'max':1000,'Nbins':10},log=True)

        See Also
        --------
        plot.h2
        """
        copy = self._ps.copy()
        a = plt.gca()

        T = np.linspace(t["min"],t["max"],t["Nbins"])

        for tt in T:
            copy.data.propagate(t=tt,verbose=True)
            plt.clf()
            a=copy.plot.h2(axis1,axis2,log=log,polar=polar,**kargs)
            a.set_title('t = %.3E fs'%tt)
            plt.pause(pause)

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
        labels=self.get_labels([axis1,axis2],kargs.get('weight','w'),kargs.get('normed',True))

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

        plt.show()

        return a

    def s2(self,axis1,axis2,weight='w',snorm=1.,log=False,polar=False,select=None):
        """
        Plot the 2d scattering plot of given axes.

        Parameters
        ----------
        axis1,axis2 : str
            name of the axes to plot
        weight : str, optional
            weight to plot. Default is w
        snorm : float, optional
            dots size normalization. Default is 1
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
        labels=self.get_labels([axis1,axis2],weight,normed=False)

        d = self._ps.data

        axis1   = d.get_axis(axis1,select=select)
        axis2   = d.get_axis(axis2,select=select)
        w       = d.get_axis(weight,select=select)

        if polar: axis1=np.radians(axis1)

        if log:
            from matplotlib.colors import LogNorm
            norm=LogNorm()
        else:
            norm=None

        s = snorm * w/max(w)

        a2 = a.scatter(axis1,axis2,c=w,s=s,norm=norm,cmap=self.cmap)

        if polar:
            pass
        else:
            a.set_xlabel(labels[0])
            a.set_ylabel(labels[1])

        plt.colorbar(a2,label=labels[2])

        plt.show()

        return a

    def h2h1(self,axis1,axis2,log=False,**kargs):
        """
        TODO
        """
        # # https://matplotlib.org/examples/pylab_examples/scatter_hist.html
        # # https://matplotlib.org/examples/axes_grid/demo_edge_colorbar.html
        # kargs1={'X':kargs.get('X',None),
        #         'erange':kargs.get('erange',[None,None]),
        #         'bwidth':kargs.get('bwidth1',None),
        #         'brange':kargs.get('brange1',[None,None])
        #         }
        #
        # kargs2={'X':kargs.get('X',None),
        #         'erange':kargs.get('erange',[None,None]),
        #         'bwidth':kargs.get('bwidth2',None),
        #         'brange':kargs.get('brange2',[None,None])
        #         }
        #
        # tmp = bool(self.autoclear)
        # if self.autoclear : self.clear()
        # self.autoclear=False
        #
        # plt.subplots_adjust(hspace=0.15,wspace=0.15)
        #
        # ax1=plt.subplot(221)
        # self.h1(axis2,log=False,reverse=True)
        # if log:ax1.set_xscale('log')
        #
        # ax2=plt.subplot(224)
        #
        # self.h1(axis1,log=log)
        #
        # ax3=plt.subplot(222,sharex=ax2,sharey=ax1)
        # self.h2(axis1,axis2,log=log,**kargs)
        # plt.setp(ax3.get_yticklabels(), visible=False)
        # plt.setp(ax3.get_xticklabels(), visible=False)
        #
        # self.autoclear=tmp
        #
        # plt.show()
        pass

    def s2h1(self,axis1,axis2,log=False):
        """
        TODO
        """
        pass

    def h3(self,axis1,axis2,axis3,s=5,wmin=0,log=False,**kargs):
        """
        Plot the 3d histogram of given axes.

        Parameters
        ----------
        axis1,axis2,axis3 : str
            Name of the axes to plot
        s : float, optional
            square sizes. Default is 5
        wmin : float, optional
            minimum weight to plot. Default is 0
        log : bool, optional
            log color scale. Default is False
        kargs : dict, optional
            Dictionnary to pass to the hist.h3 method

        See Also
        --------
        hist.h3
        """
        from mpl_toolkits.mplot3d import Axes3D
        # https://codereview.stackexchange.com/questions/62180/colorbar-for-matplotlib-3d-patch-plot
        from matplotlib import cm

        a = plt.subplot(projection='3d')

        labels=self.get_labels([axis1,axis2,axis3],kargs.get('weight','w'),kargs.get('normed',True))

        b1,b2,b3,h=self._ps.hist.h3(axis1,axis2,axis3,**kargs)

        xs,ys,zs,w = [],[],[],[]

        for i1,e1 in enumerate(h):
            for i2,e2 in enumerate(e1):
                for i3,e3 in enumerate(e2):
                    if e3>wmin:
                        xs.append(b1[i1])
                        ys.append(b2[i2])
                        zs.append(b3[i3])
                        w.append(h[i1][i2][i3])

        if log:
            c = np.log10(w)
        else:
            c = w

        a.scatter(xs=xs,ys=ys,zs=zs,c=c,s=s,marker='s')

        a.set_xlim(min(xs),max(xs))
        a.set_ylim(min(ys),max(ys))
        a.set_zlim(min(zs),max(zs))
        a.set_xlabel(labels[0],labelpad=15)
        a.set_ylabel(labels[1],labelpad=15)
        a.set_zlabel(labels[2],labelpad=15)

        m = cm.ScalarMappable(cmap=self.cmap)
        m.set_array([min(w),max(w)])

        if log:
            from matplotlib.colors import LogNorm
            norm=LogNorm()
        else:
            norm=None

        m.set_norm(norm)

        f=plt.gcf()

        f.colorbar(m,label=labels[3])

        plt.show()
        return a


    def s3(self,axis1,axis2,axis3,weight="w",snorm=1.0,log=False,select=None):
        """
        Plot the 3d histogram of given axes.

        Parameters
        ----------
        axis1,axis2,axis3 : str
            Name of the axes to plot
        weight : str, optional
            weight to plot. Default is w
        snorm : float, optional
            dots size normalization. Default is 1
        log : bool, optional
            log color scale. Default is False
        select : dict, optional
            filtering dictionnary
        """
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm

        a = plt.subplot(projection='3d')

        labels=self.get_labels([axis1,axis2,axis3],weight,normed=False)

        b1  = self._ps.data.get_axis(axis1,select=select)
        b2  = self._ps.data.get_axis(axis2,select=select)
        b3  = self._ps.data.get_axis(axis3,select=select)
        w   = self._ps.data.get_axis(weight,select=select)

        if log:
            c = np.log10(w)
        else:
            c = w
        s = snorm * w/max(w)

        a.scatter(xs=b1,ys=b2,zs=b3,c=c,s=s,marker='o')

        a.set_xlim(min(b1),max(b1))
        a.set_ylim(min(b2),max(b2))
        a.set_zlim(min(b3),max(b3))
        a.set_xlabel(labels[0],labelpad=15)
        a.set_ylabel(labels[1],labelpad=15)
        a.set_zlabel(labels[2],labelpad=15)

        m = cm.ScalarMappable(cmap=self.cmap)
        m.set_array([min(w),max(w)])

        if log:
            from matplotlib.colors import LogNorm
            norm=LogNorm()
        else:
            norm=None

        m.set_norm(norm)

        f=plt.gcf()

        f.colorbar(m,label=labels[3])

        plt.show()
        return a
