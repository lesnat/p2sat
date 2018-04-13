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
from PhaseSpace import PhaseSpaceGeneric,PhaseSpaceSmilei,PhaseSpaceGeant4


def testing(name="all"):
  from ._tests import run
  if name=="all":
    run(name="raw")
    run(name="hist")
    run(name="plot")
    run(name="heritage")
  else:
    run(name=name)
