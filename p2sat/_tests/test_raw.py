#coding:utf8
import numpy as np

from fixtures import *

def test_update():
  w,x,y,z,px,py,pz=param()
  ps=set_ps(w,x,y,z,px,py,pz)
  
  assert np.allclose(ps.raw.w,w)
  
  assert np.allclose(ps.raw.x,x)
  assert np.allclose(ps.raw.y,y)
  assert np.allclose(ps.raw.z,z)
  
  assert np.allclose(ps.raw.px,px)
  assert np.allclose(ps.raw.py,py)
  assert np.allclose(ps.raw.pz,pz)
  
def test_data():
  w,x,y,z,px,py,pz=param(size=100)
  ps=set_ps(w,x,y,z,px,py,pz)
  
  fname="fexport.p2sat"
  ps.raw.export_data(fname)
  #              w   x   y   z  px  py  pz
  ps.raw.update([0],[0],[0],[0],[0],[0],[0])
  ps.raw.import_data(fname)
  
  assert np.allclose(ps.raw.w,w)
  
  assert np.allclose(ps.raw.x,x)
  assert np.allclose(ps.raw.y,y)
  assert np.allclose(ps.raw.z,z)
  
  assert np.allclose(ps.raw.px,px)
  assert np.allclose(ps.raw.py,py)
  assert np.allclose(ps.raw.pz,pz)
  
def test_select():
  fpp=1e-7
  w,x,y,z,px,py,pz=param(size=10)
  x=np.array([5,5,5,8,7,7,1,3,8,5])
  ps=set_ps(w,x,y,z,px,py,pz)

  sw1=w[x==5]
  sw2=ps.raw.select(w,x,5,fpp=fpp)
  sw3=ps.raw.select(w,x,5.0,fpp=fpp)
  sw4=ps.raw.select('w','x',[4.9,5.1],fpp=fpp)
  
  print("w = %s"%w)
  print("x = %s"%x)
  
  
  print('sw1=%s'%sw1)
  print('sw2=%s'%sw2)
  print('sw3=%s'%sw3)
  print('sw4=%s'%sw4)
  
  assert np.allclose(sw1,sw2,rtol=fpp)
  assert np.allclose(sw2,sw3,rtol=fpp)
  assert np.allclose(sw3,sw4,rtol=fpp)

