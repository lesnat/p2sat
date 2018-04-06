#coding:utf8
import numpy as np
import matplotlib.pyplot as plt
"""
Structure :
_PhaseSpace : base object
PhaseSpaceSmilei
PhaseSpaceGeant4
PhaseSpaceTriLens
PhaseSpaceCAIN
"""
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
      """
      def update(self,w,x,y,z,px,py,pz):
        self.w  = np.array(w)
        self.x  = np.array(x)
        self.y  = np.array(y)
        self.z  = np.array(z)
        self.px = np.array(px)
        self.py = np.array(py)
        self.pz = np.array(pz)
        # 
        self.r      = np.sqrt(self.y**2+self.z**2)
        self.p      = np.sqrt(self.px**2+self.py**2+self.pz**2)
        self.theta  = np.arctan(self.py/self.px)
        self.phi    = np.arctan(self.pz/self.px) # Geometrical effect ? change -> 0 pi
        self.ekin   = (np.sqrt((self.p/0.511)**2 + 1) - 1) * 0.511
        print("Update done")

      def import_data(self,file_name,verbose=True):
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

      
      def export_data(self,file_name,title="",unit={'length':'um','momentum':'MeV/c'}):   
        print("Exporting data ...") 
        with open(file_name,'w') as f:
          f.write("# Title : %s\n"%title)
          f.write("# weight\t x ({0})\t\t y ({0})\t\t z ({0})\t\t px ({1})\t py ({1})\t pz ({1})\n".format(unit['length'],unit['momentum']))
          for i in range(len(self.w)):
            for e in [self.w[i],self.x[i],self.y[i],self.z[i],self.px[i],self.py[i],self.pz[i]]:
                f.write("{: .6E}\t".format(e))
            f.write("\n")
        print('Data succesfully exported')

    class _Hist(object):
      """
      Histograms
      """
      def __init__(self,PhaseSpace):
        self._ps=PhaseSpace
        self._r=self._ps.raw
        
      def h1(self,axis,X="all",bwidth=None,brange=[None,None]):
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
        else:
          axis= axis[self._r.x==X]
          w   = self._r.w[self._r.x==X]
        
        if bwidth is None:bwidth=axis.std()
        if brange[0] is None:brange[0]=min(axis)
        if brange[1] is None:brange[1]=max(axis)
        
        bins=np.arange(brange[0],brange[1]+bwidth,bwidth)
        
        h,b=np.histogram(axis,weights=w/bwidth,bins=bins)
        
        # Verifier b[:-1]
        return b[:-1],h
        
      def h2(self,axis1,axis2,X="all",bwidth1=None,bwidth2=None,brange1=[None,None],brange2=[None,None]):
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
        
        if X=="all":
          w   = self._r.w
        else:
          axis= axis[self._r.x==X]
          w   = self._r.w[self._r.x==X]
        
        if bwidth1 is None:bwidth1=axis1.std()
        if brange1[0] is None:brange1[0]=min(axis1)
        if brange1[1] is None:brange1[1]=max(axis1)
        bins1=np.arange(brange1[0],brange1[1]+bwidth1,bwidth1)
        if bwidth2 is None:bwidth1=axis2.std()
        if brange2[0] is None:brange2[0]=min(axis2)
        if brange2[1] is None:brange2[1]=max(axis2)
        bins2=np.arange(brange2[0],brange2[1]+bwidth2,bwidth2)
        
        h,b1,b2=np.histogram(axis1,axis2,weights=w/(bwidth1*bwidth2),bins=[bins1,bins2])
        
        # Verifier b[:-1]
        return b1[:-1],b2[:-1],h

    class _Plot(object):
      """
      Plots
      """
      def __init__(self,PhaseSpace):
        self._ps=PhaseSpace
        self._r=self._ps.raw
        self._h=self._ps.hist
        self.autoclear = False
        
      def h1(self,axis,X="all",label=["","Weight"],bins=100,log=True):
        if self.autoclear : plt.clf()
        if X=="all":
          w   = self._r.w
        else:
          axis= axis[self._r.x==X]
          w   = self._r.w[self._r.x==X]
        h,b=np.histogram(axis,weights=w,bins=bins)
        plt.step(b[:-1],h,label=label[1],where='post') # Verif
        plt.xlabel(label[0])
        if log:plt.yscale('log')
        plt.legend()

      def h2(self,axis1,axis2,X="all",label=["","","Weight"],bins=[100,100],log=False):
        if self.autoclear : plt.clf()
        if X=="all":
          w   = self._r.w
        else:
          axis1 = axis1[self._r.x==X]
          axis2 = axis2[self._r.x==X]
          w     = self._r.w[self._r.x==X]
        if log:
          from matplotlib.colors import LogNorm
          norm=LogNorm()
        else:
          norm=None
        plt.hist2d(axis1,axis2,weights=w,bins=bins,label=label[2],norm=norm)
        plt.xlabel(label[0])
        plt.ylabel(label[1])
        #plt.colorbar()
        plt.legend()
        
        
      def h1h2(self,axis1,axis2,X="all",label=["","","Weight"],bins=[100,100],log=False):
        # https://matplotlib.org/examples/pylab_examples/scatter_hist.html
        # https://matplotlib.org/examples/axes_grid/demo_edge_colorbar.html
        tmp = bool(self.autoclear)
        if self.autoclear : plt.clf()
        self.autoclear=False
        plt.subplots_adjust(hspace=0.15,wspace=0.15)
        ax1=plt.subplot(221)
        if X=="all":
            w   = self._r.w
            axis=axis1
        else:
            axis= axis1[self._r.x==X]
            w   = self._r.w[self._r.x==X]
        h,b=np.histogram(axis,weights=w,bins=bins[0])
        plt.step(h,b[:-1],label=label[2],where='post') # Verif
        plt.ylabel(label[0])
        if log:plt.xscale('log')
        plt.legend()
        ax2=plt.subplot(224)
        self.h1(axis2,X=X,label=[label[1],label[2]],bins=bins[1],log=log)
        ax3=plt.subplot(222,sharex=ax2,sharey=ax1)
        self.h2(axis1,axis2,X=X,label=['','',''],bins=bins,log=log)
        plt.setp(ax3.get_yticklabels(), visible=False)
        plt.setp(ax3.get_xticklabels(), visible=False)
        
        self.autoclear=tmp
        
      def ekin(self,X="all",bwidth=0.1): # Eventually bins=len(px)
        ekin=self._r.ekin
        bins=np.arange(min(ekin),max(ekin),bwidth)
        #self.h1(ekin,X=X,label=['ekin (MeV)','X = %s um'%X],bins=bins,log=True)
        if X=="all":
            w   = self._r.w
        else:
            ekin= ekin[self._r.x==X]
            w   = self._r.w[self._r.x==X]
        h,b=np.histogram(ekin,weights=w*1/bwidth,bins=bins)
        #h,b=np.histogram(ekin,weights=w,bins=bins)
        plt.step(b[:-1],h,label='X = %s um'%X,where='post') # Verif
        plt.xlabel('ekin (MeV)')
        plt.yscale('log')
        plt.legend()
        plt.ylabel('Number per MeV')
        
      def theta(self,X="all",bwidth=1.0):
        theta=np.degrees(self._r.theta)
        bins=np.arange(min(theta),max(theta),bwidth)
        self.h1(theta,X=X,label=['theta (deg)','X = %s um'%X],bins=bins,log=True)
        plt.ylabel('Number per bin')
      
      def phi(self,X="all",bwidth=1.0):
        phi=np.degrees(self._r.phi)
        bins=np.arange(min(phi),max(phi),bwidth)
        self.h1(phi,X=X,label=['phi (deg)','X = %s um'%X],bins=bins,log=True)
        plt.ylabel('Number per bin')

      def angle(self,X="all",bwidth=1.0):
        pass
        
      def yz(self,X="all",bwidth=10.0):
        plt.title("Number per bin")
        y = self._r.y
        z = self._r.z
        by=np.arange(-y.std(),y.std(),bwidth)
        bz=np.arange(-z.std(),z.std(),bwidth)
        #self.h1h2(y,z,X=X,label=['y (um)','z (um)','X = %s um'%X],bins=[by,bz],log=True)
        self.h1h2(y,z,X=X,label=['y (um)','z (um)',''],bins=[by,bz],log=True)
        
      def xekin(self,bwidth=0.1):
        plt.title("Number per bin")
        x = self._r.x
        ekin = self._r.ekin
        bx = [0.0]
        for e in x:
          if e not in bx:
            bx.append(e)
        bx.sort()
        bx=np.array(bx)+1e-3

        bekin=np.arange(ekin.min(),ekin.max(),bwidth)
        self.h2(x,ekin,X='all',label=['x (um)','ekin (MeV)','bin width = %s MeV'%bwidth],bins=[bx,bekin],log=True)
        
      def rekin(self,X="all",bwidth=[10.0,0.1]):
        plt.title("Number per bin")
        r = self._r.r
        ekin = self._r.ekin
        br=np.arange(0,1000,bwidth[0])
        bekin=np.arange(ekin.min(),ekin.max(),bwidth[1])
        self.h2(r,ekin,X=X,label=['r (um)','ekin (MeV)',''],bins=[br,bekin],log=True)
        
      def xn0(self):
        x = self._r.x
        bx = []
        for e in x:
          if e not in bx:
            bx.append(e)
        bx.sort()
        bx=np.array(bx)
        
        self.h1(x,label=['x (um)','Total number'],bins=bx)

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
    
    self.update(w,x,y,z,px,py,pz)
      
