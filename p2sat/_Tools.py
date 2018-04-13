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

