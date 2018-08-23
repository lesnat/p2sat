#coding:utf8
import numpy as np
import matplotlib.pyplot as plt

class _Tools(object):
  """
  Tools

  Containing various tools, such as IO into a p2sat file or fits
  """
  def __init__(self,PhaseSpace):
    self._ps=PhaseSpace

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
    t         = []
    with open(file_name,'r') as f:
      for line in f.readlines():
        if line[0]!="#":
          W,X,Y,Z,Px,Py,Pz,T=line.split()
          w.append(float(W))
          x.append(float(X))     ; y.append(float(Y))   ; z.append(float(Z))
          px.append(float(Px))   ; py.append(float(Py)) ; pz.append(float(Pz))
          t.append(float(T))
    self._ps.raw.update(w,x,y,z,px,py,pz,t,verbose=verbose)
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
      w x y z px py pz t
      w x y z px py pz t
      . . . . .  .  .  .
      . . . . .  .  .  .
      . . . . .  .  .  .
    with 7 digits precision in scientific notation

    Some text can be written if the first character of the line is a "#".

    See Also
    --------
    import_data

    TODO: add parameter 'header=True' to use header or not ?
    """
    if verbose: print("Exporting data ...")

    r=self._ps.raw

    # Opening the output file
    with open(file_name,'w') as f:

      # Write title
      f.write("# Title : %s\n"%title)

      # Write legend
      f.write("# ")
      for legend in ('weight','x (um)','y (um)','z (um)','px (MeV/c)','py (MeV/c)','pz (MeV/c)','t (fs)'):
        f.write("%-16s"%legend) # the chain is placed under 16 characters
      f.write("\n")

      # Write data
      for i in range(len(r.w)):
        for e in [r.w[i],r.x[i],r.y[i],r.z[i],r.px[i],r.py[i],r.pz[i],r.t[i]]:
            tmp="% .7E"%e # 7 digits precision with E notation
            f.write("%-16s"%tmp) # the chain is placed under 16 characters
        f.write("\n")

    if verbose: print('Data succesfully exported')

  def discretize_1(self,update=False,verbose=True,**kargs):
    """
    Discretize the particles phase space into a given number of bins.

    Parameters
    ----------
    ...


    Notes
    -----
    ...

    Examples
    --------
    ...
    """
    hn=self._ps.hist.hn
    r=self._ps.raw

    """
    if not nbins:nbins=[10]*6
    brange=[[min(r.x),max(r.x)],
            [min(r.y),max(r.y)],
            [min(r.z),max(r.z)],
            [min(r.px),max(r.px)],
            [min(r.py),max(r.py)],
            [min(r.pz),max(r.pz)]]
    bwidth=[]
    bins=[]
    for i,n in enumerate(nbins):
      bwidth.append((brange[i][1]-brange[i][0])/n)
    print(bwidth)
    """
    if verbose : print('Data discretization ...')
    #bi,hi=hn(['x','y','z','px','py','pz'],bwidth=bwidth,brange=brange,wnorm=[1.0]*6,select=select)
    bi,hi=hn([r.x,r.y,r.z,r.px,r.py,r.pz],**kargs)
    bx,by,bz,bpx,bpy,bpz=bi
    """
    hi,b=np.histogramdd([r.x,r.y,r.z,r.px,r.py,r.pz],weights=r.w,bins=nbins)

    bx=b[0]
    by=b[1]
    bz=b[2]
    bpx=b[3]
    bpy=b[4]
    bpz=b[5]
    """
    if verbose : print('Done !')
    w       = []
    x,y,z   = [],[],[]
    px,py,pz= [],[],[]

    if verbose : print('Looping over all the configurations ...')
    for ix,ex in enumerate(hi):
      if verbose: print('Step %s of %s ...'%(ix+1,len(hi)))
      for iy,ey in enumerate(ex):
        for iz,ez in enumerate(ey):
          for ipx,epx in enumerate(ez):
            for ipy,epy in enumerate(epx):
              for ipz,epz in enumerate(epy):
                if epz!=0:
                  w.append(epz)
                  x.append(bx[ix])
                  y.append(by[iy])
                  z.append(bz[iz])
                  px.append(bpx[ipx])
                  py.append(bpy[ipy])
                  pz.append(bpz[ipz])

    if verbose : print('Done !')

    if update:
      r.update(w,x,y,z,px,py,pz,verbose)
    else:
      return np.array(w),np.array(x),np.array(y),np.array(z),np.array(px),np.array(py),np.array(pz)

  def discretize(self,update=False,verbose=True,with_time=True,**kargs):
    """
    Discretize the particles phase space in a 6 or 7 D histogram.

    Parameters
    ----------
    update : bool, optional
      True for updating raw data with new discretized values. Default is False
    verbose : bool, optional
      verbosity. Default is True
    kargs
      optional keyword arguments to pass to the hist.hn function

    Notes
    -----
    This method can be used to significantly reduce disk space usage
    when saving data into output file.

    Examples
    --------
    ...

    See Also
    --------
    hist.hn
    """
    hn=self._ps.hist.hn
    r=self._ps.raw

    if verbose : print('Data discretization ...')
    #bi,hi=hn(['x','y','z','px','py','pz'],bwidth=bwidth,brange=brange,wnorm=[1.0]*6,select=select)
    if with_time:
      bi,hi=hn([r.x,r.y,r.z,r.px,r.py,r.pz,r.t],wnorm=[1.0]*7,**kargs)
      bx,by,bz,bpx,bpy,bpz,bt=bi
    else:
      bi,hi=hn([r.x,r.y,r.z,r.px,r.py,r.pz],wnorm=[1.0]*6,**kargs)
      bx,by,bz,bpx,bpy,bpz=bi
    if verbose : print('Done !')
    w       = []
    x,y,z   = [],[],[]
    px,py,pz= [],[],[]
    t       = []

    if verbose : print('Getting new configurations ...')
    hb=hi.nonzero()

    print(bx,by,bz,bpx,bpy,bpz)

    w   = hi[hb]
    x   = bx[hb[0]]
    y   = by[hb[1]]
    z   = bz[hb[2]]
    px  = bpx[hb[3]]
    py  = bpy[hb[4]]
    pz  = bpz[hb[5]]
    if with_time:
      t   = bt[hb[6]]
    else:
      t   = [0.0]*len(w)

    if verbose : print('Done !')

    if update:
      r.update(w,x,y,z,px,py,pz,t,verbose=verbose)
    else:
      return np.array(w),np.array(x),np.array(y),np.array(z),np.array(px),np.array(py),np.array(pz)

  def discretize_MT(self):
    pass

  def transformate(self,translation=(0.0,0.0,0.0),rotation=(0.0,0.0)):
    """
    Transformate the particle phase space with given translation and rotation.

    Parameters
    ----------
    ...


    Notes
    -----
    ...

    Examples
    --------
    ...
    """
    pass

  def export_TrILEns(self,path,with_time=True,verbose=True):
    """
    Export particle phase space in a TrILEns input file.

    Parameters
    ----------
    path : str
      path to the output folder
    verbose : bool, optional
      verbosity of the function. If True, a message is displayed when the data is exported
    """
    if verbose: print("Exporting data ...")

    r=self._ps.raw

    # Opening the output file
    with open(path+'prop_ph.t','w') as f:
      if with_time:
        f.write('9 1.\n')
        f.write(' poi  phx  phy  phz  pdx  pdy  pdz  gph  tim\n')
      else:
        f.write('8 1.\n')
        f.write(' poi  phx  phy  phz  pdx  pdy  pdz  gph\n')

      # Write data
      for i in range(len(r.w)):
        if with_time:
          data=[r.w[i],r.x[i],r.y[i],r.z[i],r.px[i],r.py[i],r.pz[i],r.gamma[i],r.t[i]]
        else:
          data=[r.w[i],r.x[i],r.y[i],r.z[i],r.px[i],r.py[i],r.pz[i],r.gamma[i]]

        for e in data:
            tmp="% .7E"%e # 7 digits precision with E notation
            f.write("%-16s"%tmp) # the chain is placed under 16 characters
        f.write("\n")

    if verbose: print('Data succesfully exported')


  def fit(self,axis,func_name,plot=False,polar=False,**kargs):
      """
      """
      x,w = self._ps.hist.h1(axis,**kargs)

      print(x)
      from scipy.optimize import curve_fit

      if func_name=="exp":
          f = lambda x,A,T: A*np.exp(-x/T)/T
          p0 = [sum(w),1]
      elif func_name=="gauss":
        #   f = lambda x,A,sigma: A/(np.sqrt(2*np.pi) * sigma) * np.exp(-(x)**2/(2*sigma**2))
        #   p0 = [sum(w),x.std()]
          f = lambda x,A,sigma,mu: A/(np.sqrt(2*np.pi) * sigma) * np.exp(-(x-mu)**2/(2*sigma**2))
          p0 = [sum(w),x.std(),0]
      else:
          raise NameError("Unknown func_name.")

      popt,pcov = curve_fit(f,x[:-1],w,p0=p0)

      diff = (popt[0]-sum(w))/sum(w) * 100
      print('Error on number of particles for \"{}\" fit : {:.2F} %'.format(func_name,diff))

      if plot:
          if polar:
            a=plt.gca(polar=True)
            a.plot(np.radians(x),f(x,*popt))
          else:
            a=plt.gca()
            a.plot(x,f(x,*popt))

      return x,popt
