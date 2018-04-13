#coding:utf8
import numpy as np

from fixtures import *


def test_IO_data():
  w,x,y,z,px,py,pz=param(size=100)
  ps=set_ps(w,x,y,z,px,py,pz)
  
  fname="fexport.p2sat"
  ps.tools.export_data(fname)
  #              w   x   y   z  px  py  pz
  ps.raw.update([0],[0],[0],[0],[0],[0],[0])
  ps.tools.import_data(fname)
  
  assert np.allclose(ps.raw.w,w)
  
  assert np.allclose(ps.raw.x,x)
  assert np.allclose(ps.raw.y,y)
  assert np.allclose(ps.raw.z,z)
  
  assert np.allclose(ps.raw.px,px)
  assert np.allclose(ps.raw.py,py)
  assert np.allclose(ps.raw.pz,pz)
  
