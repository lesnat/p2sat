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
    
    w,x,y,z,px,py,pz,ekin,theta,phi are instance attributes.
    
    load/import with options for different context of _PS mother class ?
    """
    def __init__(self):
        self.raw=_Raw()
        self.hist=_Hist(self)
        self.plot=_Plot(self)
      
    class _Raw(object):
      """
      Raw data
      
      Units ...
      """
      def update(self,w,x,y,z,px,py,pz):
        """
        Update class attributes
        """
        # Save values into np.array objects
        self.w  = np.array(w)
        self.x  = np.array(x)
        self.y  = np.array(y)
        self.z  = np.array(z)
        self.px = np.array(px)
        self.py = np.array(py)
        self.pz = np.array(pz)
        
        # Calculate other things ...
        self.r      = np.sqrt(self.y**2+self.z**2)
        self.p      = np.sqrt(self.px**2+self.py**2+self.pz**2)
        self.theta  = np.degrees(np.arctan(self.py/self.px))
        self.phi    = np.degrees(np.arctan(self.pz/self.px)) # Geometrical effect ? change -> 0 pi
        self.ekin   = (np.sqrt((self.p/0.511)**2 + 1) - 1) * 0.511
        
        print("Update done")

      def import_data(self,file_name,verbose=True):
        """
        Import phase space data
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
        self.update(w,x,y,z,px,py,pz)
        if verbose: print('Data succesfully imported')

      
      def export_data(self,file_name,title="",unit={'length':'um','momentum':'MeV/c'},verbose=True):
        """
        Export phase space data
        """
        if verbose: print("Exporting data ...")
        
        # Opening the output file
        with open(file_name,'w') as f:
          
          # Write title
          f.write("# Title : %s\n"%title)
          
          # Write legend
          f.write("# ")
          for legend in ('weight','x (%s)'%unit['length'],'y (%s)'%unit['length'],'z (%s)'%unit['length'],
                         'px (%s)'%unit['momentum'],'py (%s)'%unit['momentum'],'pz (%s)'%unit['momentum']):
            f.write("%-16s"%legend) # the chain is placed under 16 characters
          f.write("\n")
          
          # Write data
          for i in range(len(self.w)):
            for e in [self.w[i],self.x[i],self.y[i],self.z[i],self.px[i],self.py[i],self.pz[i]]:
                tmp="% .8E"%e # 8 digits precision with E notation
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
        
      def hn(self,axis,X="all",erange=None,bwidth=None,brange=None):
        """
        n dimension histogram
        """
        for i,ax in enumerate(axis):
          if type(ax) is str:axis[i] = eval("self._r.%s"%axis)
        
        
        for i,ax in enumerate(axis):
          if X=="all":
            w   = self._r.w
            ekin= self._r.ekin
          else:
            axis[i]= ax[self._r.x==X]
            w   = self._r.w[self._r.x==X]
            ekin= self._r.ekin[self._r.x==X]
        
        if not erange:erange=[min(ekin),max(ekin)]
        eselect=np.array([ek>erange[0] and ek<erange[1] for ek in ekin])
        w = w[eselect]
        for i,ax in enumerate(axis):
          axis=ax[eselect]
        
        if not brange:
          brange=[]
          for i,ax in enumerate(axis):
            brange.append([min(ax),max(ax)]
        
        if not bwidth:
          bwidth=[]
          for i,ax in enumerate(axis):
          bwidth.append((brange[i][1]-brange[i][0])/10.)
        
        bins=[]
        for i,ax in enumerate(axis):
          bins.append(np.arange(brange[0],brange[1]+bwidth,bwidth))
        
        h,b=np.histogramdd(axis,weights=w/np.product(bwidth),bins=bins)
        
        return *b,h
        
      def h1(self,axis,X="all",erange=[None,None],bwidth=None,brange=[None,None]):
        """
        Create and return the histo of axis
        
        Parameters
        ----------
        axis : str or np.array
          axis to hist
        X : str or float, optional
          X position on which do the hist. If "all", integrate on all the positions
        bins : int or np.array, optional
          number of bins (int) or bins on which do the hist. If None, bins=np.arange(min(axis),max(axis),axis.std())
        
        Returns
        -------
        b : np.array
          bins
        h : np.array
          histogram
        """
        if type(axis) is str:axis = eval("self._r.%s"%axis)
        
        if X=="all":
          w   = self._r.w
          ekin= self._r.ekin
        else:
          axis= axis[self._r.x==X]
          w   = self._r.w[self._r.x==X]
          ekin= self._r.ekin[self._r.x==X]
        
        if erange[0] is None:erange[0]=min(ekin)
        if erange[1] is None:erange[1]=max(ekin)
        
        eselect=np.array([ek>erange[0] and ek<erange[1] for ek in ekin])
        axis=axis[eselect]
        w = w[eselect]

        if brange[0] is None:brange[0]=min(axis)
        if brange[1] is None:brange[1]=max(axis)
        if bwidth is None:bwidth=(brange[1]-brange[0])/10.
        
        bins=np.arange(brange[0],brange[1]+bwidth,bwidth)
        
        h,b=np.histogram(axis,weights=w/bwidth,bins=bins)
        
        # Verifier b[:-1]
        h=list(h)
        h.append(0.0)
        h=np.array(h)
        return b,h
        
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
        
        print(min(axis1),max(axis1))
        print(min(axis2),max(axis2))
        
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
        
        print(bins1,bins2)
        
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
        
      def h1(self,axis,label=["",""],log=True,**kargs):
        """
        
        """
        if self.autoclear : plt.clf()
        b,h=self._h.h1(axis,**kargs)
        plt.step(b,h,'.',label=label[1],where='post') # Verif
        plt.xlim(xmin=min(b),xmax=max(b))
        plt.xlabel(label[0])
        if log:plt.yscale('log')
        plt.legend()

      def h2(self,axis1,axis2,label=["","",""],log=False,contour=True,**kargs):
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
        
      def contour(self,axis1,axis2,label=["","",""],log=False,gfilter=0.0,**kargs):
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
        
      
      def h1h2(self,axis1,axis2,label=["","",""],log=False,**kargs):
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
      
      def scatter(self,axis1,axis2,label=["","",""],log=False,**kargs):
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

      def h3(self,axis1,axis2,axis3,snorm=1.0,hmin=0.0,**kargs):
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
    self.update(w,x,y,z,px,py,pz)
    
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
    self.update(w,x,y,z,px,py,pz)


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
      
