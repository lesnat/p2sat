#coding:utf8
"""
Package Structure
=================

PhaseSpaceGeneric : mother object
PhaseSpaceXXXX : child object

Available :
PhaseSpaceSmilei
PhaseSpaceGeant4



Examples
========
Creating a random set of data

>>> import numpy as np
>>> size=1000
>>> w=np.random.uniform(low=0.0,high=100.0,size=size)
>>> x=np.random.randint(low=1.0,high=11.0,size=size)
>>> y=np.random.normal(loc=0.0,scale=1.0,size=size)
>>> z=np.random.normal(loc=0.0,scale=2.0,size=size)
>>> px=np.random.exponential(scale=1.0,size=size)
>>> py=np.random.exponential(scale=1.0,size=size)
>>> pz=np.random.exponential(scale=1.0,size=size)

Instanciate a p2sat object, let say an electron phase space, with generated data

>>> eps=PhaseSpaceGeneric()
>>> eps.raw.update(w,x,y,z,px,py,pz)

Get histograms

>>> ekin,spectrum=eps.hist.h1('ekin',bwidth=0.1,select={'x':5}) # number of e- per MeV at x==5, with a 0.1 MeV bin width

Plot some results

>>> # save bins characteristics in a dictionnary
>>> bdict = dict(bwidth1=0.5,bwidth2=1,brange1=[-5,5],brange2=[-10,10])
>>> # plots number of e- per um^2 with energy superior to 0.511 MeV, at x==5
>>> eps.plot.h2('y','z',select={'x':5,'ekin':[0.511,None]},**bdict)
>>> # add the gaussian filtered contour plot of this diag
>>> eps.plot.contour('y','z',select={'x':5,'ekin':[0.511,None]},gfilter=1.0,**bdict)






TODO
====
PhaseSpaceSmilei :
- change name & use a generic method

PhaseSpaceGeant4 :
- use a while True & try to loop over nthreads

PhaseSpaceTriLens :
- 

Code structure & names OK ?

raw :
- theta,phi schema
- theta,phi,ekin calculations

hist :
- doc hn
- use nbins OK ?

plot :
- doc

tools :?
- fit MB, gaussian, MJ
- IO(file_name,mode='r',title="")
"""

import numpy as np
import matplotlib.pyplot as plt


class PhaseSpaceGeneric(object):
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
      - ekin is defined as :math:`(\sqrt{(p/m_e c)^2+1}-1) \\times m_e c^2` *
      - theta is defined as :math:`\\arctan{p_y/p_x}`
      - phi is defined (yet) as :math:`\\arctan{p_z/p_x}`
      
      
      * detail of the calculus can be found at ... TODO
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
          filtering value/range (value if int, range if float or list/tuple). If a frange element is None, the minimum/maximum value is taken
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
        >>> w = select(w,ekin,[0.511,None]) # Select all the w with :math:`ekin \in [0.511,+\infty] MeV`
        
        If frange is a list/tuple or a float, the filtering is done with a fpp precision
        """
        if type(axis) is str: axis=eval("self.%s"%axis)
        if type(faxis) is str: faxis=eval("self.%s"%faxis)
        
        if type(frange) is list or type(frange) is tuple:
          if frange[0] is None: frange[0]=min(axis)
          if frange[1] is None: frange[1]=max(axis)
          filtr=np.array([x>frange[0]*(1-fpp) and x<frange[1]*(1+fpp) for x in faxis])
          axis=axis[filtr]
        elif type(frange) is int:
          axis=axis[faxis==frange]
        elif type(frange) is float:
          axis=self.select(axis,faxis,frange=[frange,frange])
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
          verbosity of the function. If True, a message is displayed when the data is imported
        
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
          verbosity of the function. If True, a message is displayed when the data is exported
        
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
        import_data
        
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
        it reduce bwidth to fit brange[1]-brange[0] with a int nb of bins
        
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
            for i,ax in enumerate(axis):
              axis[i]=self._r.select(ax,key,val)
        
        # Define default bin range
        if not brange: brange=[[None,None]]*len(axis) 
        for i,br in enumerate(brange):
          if br[0] is None:brange[i][0]=min(axis[i])
          if br[1] is None:brange[i][1]=max(axis[i])
        
        # Define default bin width + number of bins
        nbins=[]
        if not bwidth: bwidth=[None]*len(axis)
        for i,bw in enumerate(bwidth):
          blen=brange[i][1] - brange[i][0]
          if bw is None: bwidth[i]=blen/10.
          nbins.append(np.ceil(blen/bwidth[i] + 1))
        
        # Construct the bins list
        bins=[]
        for i,ax in enumerate(axis):
          bins.append(np.linspace(brange[i][0],brange[i][1],int(nbins[i])))
          #bins.append(np.arange(brange[i][0],brange[i][1],bwidth[i]))
        
        # Define weight normalization
        if not wnorm: wnorm=[None]*len(axis)
        for i,wn in enumerate(wnorm):
          if wn is None: wnorm[i]=bwidth[i]
        
        # Calculate the multi dimensional histo, normalized by wnorm
        h,b=np.histogramdd(axis,weights=w/np.product(wnorm),bins=bins)
        
        # Return the bins and histo
        return b,h
        
      def h1(self,axis,bwidth=None,brange=None,select=None):
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
        
        b,h=self.hn([axis],bwidth=[bwidth],brange=[brange],select=select)
        
        # Verifier b[:-1]
        h=list(h)
        h.append(0.0)
        h=np.array(h)
        return b[0],h
        
      def h2(self,axis1,axis2,bwidth1=None,bwidth2=None,brange1=None,brange2=None,select=None):
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
        
        b,h=self.hn([axis1,axis2],bwidth=[bwidth1,bwidth2],brange=[brange1,brange2],select=select)
        
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
    
class PhaseSpaceSmilei(PhaseSpaceGeneric):
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


class PhaseSpaceGeant4(PhaseSpaceGeneric):
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
      
