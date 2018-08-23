#coding:utf8
import numpy as np

class _Extract(object):
  """
  Import data from a file.

  Notes
  -----
  If you want to add a method to import data from another code, you must proceed as follow :

  - Add a method to this object, with the name of your code. It must contains the keyword `self` as a first argument (because of object-oriented paradigm), and all the other parameters you need
  - Get the data from your file and put it in lists or numpy arrays, one line describing one particle
  - Call the `update` method of `data` sub-object (access via `self._ps.data.update`)
  - Please write a documentation and share !

  You can copy-paste the `txt` method to have a basic example of file import.
  """
  def __init__(self,PhaseSpace):
    self._ps=PhaseSpace

  def txt(self,file_name,sep=",",verbose=True):
    """
    Load particle phase space from a text file.

    Parameters
    ----------
    file_name : str
      name of the input file
    sep : str
      character used to separate values. Default is ','
    verbose : bool, optional
      verbosity of the function. If True, a message is displayed when the data is imported

    See Also
    --------
    export.txt
    """
    if verbose: print("Importing data ...")
    w         = []
    x,y,z     = [],[],[]
    px,py,pz  = [],[],[]
    t         = []
    with open(file_name,'r') as f:
      for line in f.readlines():
        if line[0]!="#":
          W,X,Y,Z,Px,Py,Pz,T=line.split(sep)
          w.append(float(W))
          x.append(float(X))     ; y.append(float(Y))   ; z.append(float(Z))
          px.append(float(Px))   ; py.append(float(Py)) ; pz.append(float(Pz))
          t.append(float(T))
    self._ps.data.update(w,x,y,z,px,py,pz,t,verbose=verbose)
    if verbose: print('Data succesfully imported')

  def Smilei_Screen_1d(self,Screen,xnorm,wnorm,tnorm,X=0):
    w         = []
    x,y,z     = [],[],[]
    px,py,pz  = [],[],[]
    t         = []

    time= Screen().get()['times']

    Px  = Screen().get()['px'] * 0.511
    Py  = Screen().get()['py'] * 0.511

    wNorm = wnorm/(0.511**2)
    wNorm *= (max(Px)-min(Px))/len(Px)
    wNorm *= (max(Py)-min(Py))/len(Py)

    cdata=np.array([[0.]*len(Px)]*len(Py))

    print("Extracting screen data ...")
    for it,et in enumerate(time):
      ldata = cdata
      cdata = Screen(timesteps=et).getData()[0]
      for ipx,epx in enumerate(cdata):
        for ipy,epy in enumerate(epx):
          depy = epy-ldata[ipx][ipy]
          if depy!=0.:
            w.append(depy*wNorm)
            x.append(X)
            y.append(0.)
            z.append(0.)
            px.append(Px[ipx])
            py.append(Py[ipy])
            pz.append(0.)
            t.append(et*tnorm)


    pz = [0.0] * len(w)
    x = [X] * len(w)
    z = [0.0] * len(w)
    print("Data succesfully imported")
    self._ps.data.update(w,x,y,z,px,py,pz,t)

  def Geant4_csv(self,file_name,nthreads=1,verbose=True):
    """
    Extract simulation results from a Geant4 NTuple csv output file

    Parameters
    ----------
    file_name : str
      name of the output file. If it ends with '*_t0.*', the number '0' will be replaced by the number of the current thread
    nthreads : int
      total number of threads to consider
    verbose : bool, optional
      verbosity

    Notes
    -----
    The Geant4 NTuple format should be
      w,x,y,z,px,py,pz
      . . . . .  .  .
      . . . . .  .  .
      . . . . .  .  .

    Examples
    --------
    >>> eg = p2sat.PhaseSpaceGeant4()
    >>> eg.extract("../Geant4/testem_nt_electron_t*.csv",nthreads=10)

    TODO
    ----
    - while True + try/except to loop over nthreads ?
    """
    data = []
    fext = file_name.split('.')[-1]   # File extension
    fbase= file_name[:-(len(fext)+1)] # File base name

    for thread in range(0,nthreads):
      fname=fbase[:-1]+str(thread)+"."+fext
      if verbose:print("Extracting %s ..."%fname)

      if fext=="csv":
        with open(fname,'r') as f:
          for line in f.readlines():
            if line[0]!='#':
              for e in line.split(','):
                data.append(float(e))
      elif fext=="xml":
        from lxml import etree
        with etree.parse(fname) as tree:
          for entry in tree.xpath('/aida/tuple/rows/row/entry'):
            data.append(float(entry.get('value')))
      else:
        raise NameError("Unknown file extension : %s"%fext)

    w   = data[0::8]
    x   = data[1::8]
    y   = data[2::8]
    z   = data[3::8]
    px  = data[4::8]
    py  = data[5::8]
    pz  = data[6::8]
    t   = data[7::8]
    if verbose:print("Done !")

    self._ps.data.update(w,x,y,z,px,py,pz,t,verbose)

  def TrILEns_output(self,path,specie,verbose=True):
    """
    Extract simulation results from a TrILEns output.txt file

    Parameters
    ----------
    path : str
      simulation path
    specie : str
      specie to find in the output. The specie name must be in plural form (i.e 'electrons' or 'positrons')
    verbose : bool, optional
      verbosity

    Notes
    -----
    ...

    Examples
    --------
    >>> et = p2sat.PhaseSpaceTrILEns()
    >>> et.extract("../TrILEns/",specie="positrons")

    TODO
    ----
    - Check Interface.f90 at line ~ 120 -> condition for e-/e+ & not working with photons
    """
    if verbose:print("Extracting {} data from {}output.txt ...".format(specie,path))

    w         = []
    x,y,z     = [],[],[]
    px,py,pz  = [],[],[]
    t         = []

    is_correct_specie=False

    with open(path+'output.txt','r') as f:
      for _ in range(34):
        f.readline()

      for line in f.readlines():
        try:
          W,X,Y,Z,Px,Py,Pz,Gamma,Chi=line.split()
          if is_correct_specie:
            w.append(float(W))
            x.append(float(X))     ; y.append(float(Y))   ; z.append(float(Z))
            px.append(float(Px))   ; py.append(float(Py)) ; pz.append(float(Pz))
            t.append(0.)
        except ValueError:
          if specie in line.split():
            is_correct_specie = True
          else:
            is_correct_specie = False

    if verbose:print("Done !")

    self._ps.data.update(w,x,y,z,px,py,pz,t,verbose=verbose)
