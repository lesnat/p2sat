#coding:utf8
import numpy as np

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
    with open(file_name,'r') as f:
      for line in f.readlines():
        if line[0]!="#":
          W,X,Y,Z,Px,Py,Pz=line.split()
          w.append(float(W))
          x.append(float(X))     ; y.append(float(Y))   ; z.append(float(Z))
          px.append(float(Px))   ; py.append(float(Py)) ; pz.append(float(Pz))
    self._ps.raw.update(w,x,y,z,px,py,pz,verbose=verbose)
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

    r=self._ps.raw

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
      for i in range(len(r.w)):
        for e in [r.w[i],r.x[i],r.y[i],r.z[i],r.px[i],r.py[i],r.pz[i]]:
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

  def discretize_2(self,update=False,verbose=True,**kargs):
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
    print(r.w)
    if verbose : print('Data discretization ...')
    #bi,hi=hn(['x','y','z','px','py','pz'],bwidth=bwidth,brange=brange,wnorm=[1.0]*6,select=select)
    bi,hi=hn([r.x,r.y,r.z,r.px,r.py,r.pz],wnorm=[1.0]*6,**kargs)
    bx,by,bz,bpx,bpy,bpz=bi
    print(hi)
    if verbose : print('Done !')
    w       = []
    x,y,z   = [],[],[]
    px,py,pz= [],[],[]

    if verbose : print('Getting new configurations ...')
    hb=hi.nonzero()
    print(hb)

    w   = hi[hb]
    x   = bx[hb[0]]
    y   = by[hb[1]]
    z   = bz[hb[2]]
    px  = bpx[hb[3]]
    py  = bpy[hb[4]]
    pz  = bpz[hb[5]]

    if verbose : print('Done !')

    if update:
      r.update(w,x,y,z,px,py,pz,verbose)
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
