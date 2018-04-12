#coding:utf8
from __future__ import absolute_import
import numpy as np
from ..PhaseSpace import PhaseSpaceGeneric

__all__=['param','set_ps','get_ps']


def param(size=1000):
  w=np.random.uniform(low=0.0,high=100.0,size=size)
  
  x=np.random.randint(low=1.0,high=11.0,size=size)
  y=np.random.normal(loc=0.0,scale=1.0,size=size)
  z=np.random.normal(loc=0.0,scale=2.0,size=size)
  
  px=np.random.exponential(scale=1.0,size=size)
  py=np.random.exponential(scale=1.0,size=size)
  pz=np.random.exponential(scale=1.0,size=size)
  
  return w,x,y,z,px,py,pz
  
def set_ps(w,x,y,z,px,py,pz):
  ps=PhaseSpaceGeneric()
  ps.raw.update(w,x,y,z,px,py,pz,verbose=False)
  return ps
  
def get_ps(size=1000):
  return set_ps(*param(size=size))

