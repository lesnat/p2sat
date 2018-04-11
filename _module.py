#coding:utf8
"""
Package Structure
=================

_PhaseSpace : base object
PhaseSpaceXXXX : child object

Available :
PhaseSpaceSmilei
PhaseSpaceGeant4

TODO : 
PhaseSpaceTriLens
PhaseSpaceCAIN
"""

import numpy as np
import matplotlib.pyplot as plt


class _PhaseSpace(object):
    """
    Base class for particle phase-space analysis.
    
    Attributes
    ----------
    extract : "polymorphic" method
      load phase space from a simulation file
    raw : sub-object
      contains raw data and methods to manipulate it, such as import/export into a file or value filtering method. Appropriate place to IO ?
    hist : sub-object
      contains methods to make histograms from raw data
    plot : sub-object
      contains methods to plot histos
    
    Notes
    -----
    See sub-objects documentation for more informations
    """
    def __init__(self):
        self.raw=self._Raw()
        self.hist=self._Hist(self)
        self.plot=self._Plot(self)
      
    def extract(self):
      """
      Extract raw data from a simulation file.
      
      Notes
      -----
      This abstract method must be overwritten by the _PhaseSpace child classes
      """
      raise NotImplementedError
      
    class _Raw(object):
      """
      Class containing raw data and methods to manipulate it.
      
      Attributes
      ----------
      w : numpy.ndarray
        particle statistical weight
      x,y,z : numpy.ndarray
        particle x,y,z position in um
      r : numpy.ndarray
        absolute distance to the x axis in um
      px,py,pz : numpy.ndarray
        particle momentum in x,y,z direction in MeV/c
      p : numpy.ndarray
        absolute particle momentum in MeV/c
      ekin : numpy.ndarray
        particle energy in MeV
      theta : numpy.ndarray
        angle between px and py in degree
      phi : numpy.ndarray
        angle between ??? in degree
      
      Notes
      -----
      As all the calculations are done with the previously defined units,
      the input data might be firstly converted to those units.
      
      Calculations :
      - r is defined as :math:`\sqrt{y^2+z^2}`
      - p is defined as :math:`\sqrt{p_x^2+p_y^2+p_z^2}`
      - ekin is defined as :math:`(\sqrt{(p/m_e c)^2+1}-1) \\times m_e c^2`
      - theta is defined as :math:`\\arctan{p_y/p_x}`
      - phi is defined (yet) as :math:`\\arctan{p_z/p_x}`
      """
      def select(self,axis,faxis,frange,fpp=1e-7):
        """
        Filter an axis with a value/range on another axis.
        
        Parameters
        ----------
        axis : str or numpy.ndarray
          axis to filter
        faxis : str or numpy.ndarray
          filtering axis
        frange : int, float, list/tuple of 2 float
          filtering value/range (value if int, range if float or list/tuple)
        fpp : float, optional
          relative floating point precision. Default is 1e-7
          
        Returns
        -------
        axis : numpy.ndarray
          filtered axis
          
        Examples
        --------
        
        It is possible to filter by an int value
        
        >>> w = np.random.uniform(low=0.,high=10.,size=10)
        >>> x = np.array([1,3,3,3,7,9,5,3,7,3])
        >>> w = select(w,x,3) # Select all the w satisfying x==3
        
        or filter by a range
        
        >>> w = np.random.uniform(low=0.,high=10.,size=1000)
        >>> ekin = np.random.exponential(scale=3.0,size=1000)
        >>> w = select(w,ekin,[0.511,None]) # Select all the w with :math:`ekin \in [0.511,+\infty]`
        
        If frange is a list/tuple or a float, the filtering is done with a fpp precision
        """
        if type(axis) is str: axis=eval("self.%s"%axis)
        if type(faxis) is str: faxis=eval("self.%s"%faxis)
        
        if type(frange) is list or type(frange) is tuple:
          select=np.array([x>frange[0]-fpp and x<frange[1]+fpp for x in faxis])
          axis=axis[select]
        elif type(frange) is int:
          axis=axis[faxis==frange]
        elif type(frange) is float:
          axis=self._filter(axis,faxis,frange=(frange-fpp*frange,frange+fpp*frange))
        else:
          raise TypeError('frange type must be int/float or list/tuple of 2 float.')
        
        return axis
        
      def update(self,w,x,y,z,px,py,pz,verbose=True):
        """
        Update class attributes with new values.
        
        Parameters
        ----------
        w,x,y,z,px,py,pz : list or numpy.ndarray
          particle phase space. More information can be found in raw object documentation
        verbose : bool
          verbosity of the function. If True, a message is displayed when the attributes are loaded in memory
        """
        # Save values into np.array objects
        self.w  = np.array(w)
        self.x  = np.array(x)
        self.y  = np.array(y)
        self.z  = np.array(z)
        self.px = np.array(px)
        self.py = np.array(py)
        self.pz = np.array(pz)
        
        # Calculate other parameters from it
        self.r      = np.sqrt(self.y**2+self.z**2)
        self.p      = np.sqrt(self.px**2+self.py**2+self.pz**2)
        self.theta  = np.degrees(np.arctan(self.py/self.px))
        self.phi    = np.degrees(np.arctan(self.pz/self.px)) # Geometrical effect ? change -> 0 pi
        self.ekin   = (np.sqrt((self.p/0.511)**2 + 1) - 1) * 0.511
        
        if verbose: print("Update done")

      def import_data(self,file_name,verbose=True):
        """
        Import particle phase space from a p2sat file.
        
        Parameters
        ----------
        file_name : str
          name of the input file
        verbose : bool, optional
          verbosity
        
        See Also
        --------
        export_data
        """
        if verbose: print("Importing data ...")
        w         = []
        x,y,z     = [],[],[]
        px,py,pz  = [],[],[]
        with open(file_name,'r') as f:
          for line in f.readlines():
            if line[0]!="#":
              W,X,Y,Z,Px,Py,Pz=line.split()
              w.append(float(W))
              x.append(float(X))     ; y.append(float(Y))   ; z.append(float(Z))
              px.append(float(Px))   ; py.append(float(Py)) ; pz.append(float(Pz))
        self.update(w,x,y,z,px,py,pz,verbose=verbose)
        if verbose: print('Data succesfully imported')

      
      def export_data(self,file_name,title="",verbose=True):
        """
        Export particle phase space in a p2sat file.
        
        Parameters
        ----------
        file_name : str
          name of the output file
        title : str, optional
          title of the file
        verbose : bool, optional
          verbosity
        
        Notes
        -----
        The format in the output file is
        # title
        # legend
          w x y z px py pz
          w x y z px py pz
          . . . . .  .  .
          . . . . .  .  .
          . . . . .  .  .
        with 7 digits precision in scientific notation
        
        Some text can be written if the first character of the line is a "#".
        
        See Also
        --------
        import_p2sat
        
        TODO: add parameter 'header=True' to use header or not ?
        """
        if verbose: print("Exporting data ...")
        
        # Opening the output file
        with open(file_name,'w') as f:
          
          # Write title
          f.write("# Title : %s\n"%title)
          
          # Write legend
          f.write("# ")
          for legend in ('weight','x (um)','y (um)','z (um)','px (MeV/c)','py (MeV/c)','pz (MeV/c)'):
            f.write("%-16s"%legend) # the chain is placed under 16 characters
          f.write("\n")
          
          # Write data
          for i in range(len(self.w)):
            for e in [self.w[i],self.x[i],self.y[i],self.z[i],self.px[i],self.py[i],self.pz[i]]:
                tmp="% .7E"%e # 7 digits precision with E notation
                f.write("%-16s"%tmp) # the chain is placed under 16 characters
            f.write("\n")
        
        if verbose: print('Data succesfully exported')

    class _Hist(object):
      """
      Histograms
      """
      def __init__(self,PhaseSpace):
        self._ps=PhaseSpace
        self._r=self._ps.raw
        
      def hn(self,axis,bwidth=None,brange=None,wnorm=None,select=None):
        """
        Create and return the n-dimensional histo of axis list.
        
        Parameters
        ----------
        axis : list of str/np.array
          list of axis to hist
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
        
        Example
        -------
        
        >>> hn(['x'],bwidth=[50],brange=[[0,1000]],wnorm=[1.0],select={'ekin':(0.511,None)})
        
        returns the number of particles with :math:`ekin \in [0.511, +\infty] MeV` in function of x
        wnorm=[1.0] to not divide nb of particles by bin width (otherwise number per um)
        
        >>> hn(['r','ekin'],bwidth=[10.0,0.1],brange=[[0,1000],[0.1,50.0]],select={'x':150})
        
        returns a number of e- per um per MeV at x=150 um
        """
        # Get the axis from a str if needed
        for i,ax in enumerate(axis):
          if type(ax) is str:axis[i] = eval("self._r.%s"%ax)
        
        # Get a shortcut to particle statistical weight
        w   = self._r.w
        
        # Filter the data if needed
        if select is not None:
          for key,val in select.items():
            w = self._r.select(w,key,val)
            axis=self._r.select(axis,key,val)
        
        # Define default bin range
        if not brange: brange=[[None,None]]*len(axis) 
        for i,br in enumerate(brange):
          if br[0] is None:brange[i][0]=min(axis[i])
          if br[1] is None:brange[i][1]=max(axis[i])
        
        # Define default bin width
        if not bwidth: bwidth=[None]*len(axis)
        for i,bw in enumerate(bwidth):
          if bw is None: bwidth[i]=(brange[i][1]-brange[i][0])/10.
        
        # Calculate number of bins
        
        # Construct the bins list
        bins=[]
        for i,ax in enumerate(axis):
          #bins.append(np.linspace(brange[i][0],brange[i][1],nbins[i]))
          bins.append(np.arange(brange[i][0],brange[i][1],bwidth[i]))
        
        # Define weight normalization
        if not wnorm: wnorm=[None]*len(axis)
        for i,wn in enumerate(wnorm):
          if wn is None: wnorm[i]=bwidth[i]
        
        # Calculate the multi dimensional histo, normalized by wnorm
        h,b=np.histogramdd(axis,weights=w/np.product(wnorm),bins=bins)
        
        # Return the bins and histo
        return b,h
        
      def h1(self,axis,X=None,erange=None,bwidth=None,brange=None):
        """
        Create and return the 1 dimensional histogram of given axis.
        
        Parameters
        ----------
        axis : str or np.array
          axis to hist
        
        ...
        
        Returns
        -------
        b : np.array
          bins
        h : np.array
          histogram
          
        Notes
        -----
        TODO: append 0 at the end ?
        
        See Also
        --------
        hn,h2,h3
        """
        b,h=self.hn([axis],X=X,erange=erange,bwidth=[bwidth],brange=[brange])
        
        # Verifier b[:-1]
        h=list(h)
        h.append(0.0)
        h=np.array(h)
        return b[0],h
        
      def h2(self,axis1,axis2,X="all",erange=[None,None],bwidth1=None,bwidth2=None,brange1=[None,None],brange2=[None,None]):
        """
        Create and return the histo of axis
        
        Parameters
        ----------
        axis1,axis2 : str or np.array
          axis to hist
        X : str or float, optional
          X position on which do the hist. If "all", integrate on all the positions
        bins : int or np.array, optional
          number of bins (int) or bins on which do the hist. If None, bins=np.arange(min(axis),max(axis),axis.std())
        
        Returns
        -------
        b1,b2 : np.array
          bins
        h : np.array
          histogram
          
          
        TODO: NE PAS AFFECTER DE VALEURS AUX PARAMETRES -> Reaffecte la valeur par défaut si utilisé dans un module
        """
        if type(axis1) is str:axis1 = eval("self._r.%s"%axis1)
        if type(axis2) is str:axis2 = eval("self._r.%s"%axis2)
        
        if X=="all":
          w   = self._r.w
          ekin= self._r.ekin
        else:
          axis1= axis1[self._r.x==X]
          axis2= axis2[self._r.x==X]
          w   = self._r.w[self._r.x==X]
          ekin= self._r.ekin[self._r.x==X]
          
        if erange[0] is None:erange[0]=min(ekin)
        if erange[1] is None:erange[1]=max(ekin)
        
        eselect=np.array([ek>erange[0] and ek<erange[1] for ek in ekin])
        axis1=axis1[eselect]
        axis2=axis2[eselect]
        w = w[eselect]
        
        if brange1[0] is None:brange1[0]=min(axis1)
        if brange1[1] is None:brange1[1]=max(axis1)
        if bwidth1 is None:bwidth1=(brange1[1]-brange1[0])/10.
        bins1=np.arange(brange1[0],brange1[1]+bwidth1,bwidth1)
        if brange2[0] is None:brange2[0]=min(axis2)
        if brange2[1] is None:brange2[1]=max(axis2)
        if bwidth2 is None:bwidth2=(brange2[1]-brange2[0])/10.
        bins2=np.arange(brange2[0],brange2[1]+bwidth2,bwidth2)
        
        h,b1,b2=np.histogram2d(axis1,axis2,weights=w/(bwidth1*bwidth2),bins=[bins1,bins2])
        
        # Verifier b[:-1]
        #return b1[:-1],b2[:-1],h
        return b1,b2,h
        
      def h3(self,axis1,axis2,axis3,X="all",erange=[None,None],bwidth1=None,bwidth2=None,bwidth3=None,brange1=[None,None],brange2=[None,None],brange3=[None,None]):
        """
        Create and return the histo of axis
        
        Parameters
        ----------
        axis1,axis2 : str or np.array
          axis to hist
        X : str or float, optional
          X position on which do the hist. If "all", integrate on all the positions
        bins : int or np.array, optional
          number of bins (int) or bins on which do the hist. If None, bins=np.arange(min(axis),max(axis),axis.std())
        
        Returns
        -------
        b1,b2 : np.array
          bins
        h : np.array
          histogram
        """
        if type(axis1) is str:axis1 = eval("self._r.%s"%axis1)
        if type(axis2) is str:axis2 = eval("self._r.%s"%axis2)
        if type(axis3) is str:axis3 = eval("self._r.%s"%axis3)
        
        if X=="all":
          w   = self._r.w
          ekin= self._r.ekin
        else:
          axis1= axis1[self._r.x==X]
          axis2= axis2[self._r.x==X]
          axis3= axis3[self._r.x==X]
          w   = self._r.w[self._r.x==X]
          ekin= self._r.ekin[self._r.x==X]
          
        if erange[0] is None:erange[0]=min(ekin)
        if erange[1] is None:erange[1]=max(ekin)
        
        eselect=np.array([ek>erange[0] and ek<erange[1] for ek in ekin])
        axis1=axis1[eselect]
        axis2=axis2[eselect]
        axis3=axis3[eselect]
        w = w[eselect]
        
        if brange1[0] is None:brange1[0]=min(axis1)
        if brange1[1] is None:brange1[1]=max(axis1)
        if bwidth1 is None:bwidth1=(brange1[1]-brange1[0])/10.
        bins1=np.arange(brange1[0],brange1[1]+bwidth1,bwidth1)
        if brange2[0] is None:brange2[0]=min(axis2)
        if brange2[1] is None:brange2[1]=max(axis2)
        if bwidth2 is None:bwidth2=(brange2[1]-brange2[0])/10.
        bins2=np.arange(brange2[0],brange2[1]+bwidth2,bwidth2)
        if brange3[0] is None:brange3[0]=min(axis3)
        if brange3[1] is None:brange3[1]=max(axis3)
        if bwidth3 is None:bwidth3=(brange3[1]-brange3[0])/10.
        bins3=np.arange(brange3[0],brange3[1]+bwidth3,bwidth3)
        
        h,b=np.histogramdd([axis1,axis2,axis3],weights=w/(bwidth1*bwidth2*bwidth3),bins=[bins1,bins2,bins3])
        
        # Verifier b[:-1]
        #return b[0][:-1],b[1][:-1],b[2][:-1],h
        return b[0],b[1],b[2],h

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
        g1,g2=np.meshgrid(b1,b2)
        a1=plt.pcolormesh(g1,g2,h.T,norm=norm) # Voire .T
        
        plt.xlim(xmin=min(b1),xmax=max(b1))
        plt.ylim(ymin=min(b2),ymax=max(b2))
        
        plt.xlabel(label[0])
        plt.ylabel(label[1])
        plt.title(label[2])
        ca1=plt.colorbar(a1,orientation='horizontal')
        
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
        g1,g2=np.meshgrid(b1,b2)
        #if gfilter>0.0:
        from scipy.ndimage.filters import gaussian_filter
        a2=plt.contour(g1[:-1,:-1],g2[:-1,:-1],gaussian_filter(h.T,gfilter),
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
        g1,g2=np.meshgrid(b1,b2)
        axh2.pcolormesh(g1,g2,h.T) # Voire .T
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
        g1,g2,g3=np.meshgrid(b1,b2,b3)
        
        tmp=np.array(h)
        
        for i1,e1 in enumerate(h):
          for i2,e2 in enumerate(e1):
            for i3,e3 in enumerate(e2):
              if e3<hmin:
                tmp[i1][i2][i3]=0.0
              else:
                tmp[i1][i2][i3]=e3
        
        ax.scatter3D(g1,g2,g3,s=snorm*tmp.T,cmap='hot')

#########################################################################################
    
class PhaseSpaceSmilei(_PhaseSpace):
  def __init__(self):
    self.raw = self._Raw()
    self.hist = self._Hist(self)
    self.plot = self._Plot(self)
    
  def extract_1d(self,Screen,timestep,xnorm,wnorm=1.0,X=0):
    w         = []
    x,y,z     = [],[],[]
    px,py,pz  = [],[],[]
    
    data= Screen(timesteps=timestep).getData()[0]
    Px  = Screen().get()['px'] * 0.511
    Py  = Screen().get()['py'] * 0.511
    
    try:
      Y = Screen().get()['y'] * xnorm
    except:
      Y = [0.0]
      data = [data]
    
    wNorm = wnorm/(0.511**2)
    wNorm *= (max(Px)-min(Px))/len(Px)
    wNorm *= (max(Py)-min(Py))/len(Py)
    
    # redéfinir l'ordre dans namelist ?
    print("Extracting screen data ...")
    for iy,ey in enumerate(data):
      for ipx,epx in enumerate(ey):
        for ipy,epy in enumerate(epx):
          if epy!=0:
            y.append(Y[iy])
            px.append(Px[ipx])
            py.append(Py[ipy])
            w.append(epy*wNorm)

    pz = [0.0] * len(w)
    x = [X] * len(w)
    z = [0.0] * len(w)
    print("Data succesfully imported")
    self.raw.update(w,x,y,z,px,py,pz)
    
  def extract_2d(self,Screen,timestep,xnorm,wnorm=1.0,X=0):
    w         = []
    x,y,z     = [],[],[]
    px,py,pz  = [],[],[]
    
    data= Screen(timesteps=timestep).getData()[0]
    Px  = Screen().get()['px'] * 0.511
    Py  = Screen().get()['py'] * 0.511
    Pz  = Screen().get()['pz'] * 0.511


    
    try:
      Y = Screen().get()['y'] * xnorm
    except:
      Y = [0.0]
      data = [data]
    
    wNorm  = wnorm/(0.511**3)
    wNorm *= (max(Px)-min(Px))/len(Px)
    wNorm *= (max(Py)-min(Py))/len(Py)
    wNorm *= (max(Pz)-min(Pz))/len(Pz)
    # redéfinir l'ordre dans namelist ?
    print("Extracting screen data ...")
    for iy,ey in enumerate(data):
      for ipx,epx in enumerate(ey):
        for ipy,epy in enumerate(epx):
          for ipz,epz in enumerate(epy):
            if epz!=0:
              y.append(Y[iy])
              px.append(Px[ipx])
              py.append(Py[ipy])
              pz.append(Pz[ipz])
              w.append(epz*wNorm)
              
    x = [X] * len(w)
    z = [0.0] * len(w)
    print("Data succesfully imported")
    self.raw.update(w,x,y,z,px,py,pz)


class PhaseSpaceGeant4(_PhaseSpace):
  def __init__(self):
    self.raw = self._Raw()
    self.hist = self._Hist(self)
    self.plot = self._Plot(self)
    
  def extract(self,file_name,nthreads=1):
    data = []
    fext = file_name[-4:]
    
    if fext=="root":
      pass
    else:
      for thread in range(0,nthreads):
        fname=file_name[:-5]+"%s"%thread+fext
        print("Extracting %s ..."%fname)

        if fext==".csv":
          with open(fname,'r') as f:
            for line in f.readlines():
              if line[0]!='#':
                for e in line.split(','):
                  data.append(float(e))
        
        if fext==".xml":
          from lxml import etree      
          tree = etree.parse(fname)
          for entry in tree.xpath('/aida/tuple/rows/row/entry'):
            data.append(float(entry.get('value')))
    
    w   = data[0::7]
    x   = data[1::7]
    y   = data[2::7]
    z   = data[3::7]
    px  = data[4::7]
    py  = data[5::7]
    pz  = data[6::7]
    print("Data succesfully imported")
    
    self.raw.update(w,x,y,z,px,py,pz)
      
