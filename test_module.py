#coding:utf8
import _module
import numpy as np

def param():
  size=100000
  
  w=np.random.uniform(low=0.0,high=100.0,size=size)
  
  x=np.random.exponential(scale=1.0,size=size)
  y=np.random.normal(loc=0.0,scale=1.0,size=size)
  z=np.random.normal(loc=0.0,scale=2.0,size=size)
  
  px=np.random.exponential(scale=1.0,size=size)
  py=np.random.exponential(scale=1.0,size=size)
  pz=np.random.exponential(scale=1.0,size=size)
  
  return w,x,y,z,px,py,pz
  
def set_ps(w,x,y,z,px,py,pz):
  ps=_module._PhaseSpace()
  ps.raw.update(w,x,y,z,px,py,pz,verbose=False)
  return ps
  
def get_ps():
  return set_ps(*param())

  
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
  w,x,y,z,px,py,pz=param()
  ps=set_ps(w,x,y,z,px,py,pz)
  
  fname="test_export.dat"
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
  
  

def test_hist_hn_default():
  ps=get_ps()
  b,h=ps.hist.hn(['px'])
  b = b[0]
  
  # Test default bin length
  assert len(b) == 10
  
  # Test default brange
  assert np.isclose(b[-1],max(ps.raw.x))
  assert np.isclose(b[0],min(ps.raw.x))
  
  # Test erange selection
  
def test_hist_hn_erange():
  ps=get_ps()

def test_hist_h1():
  w,x,y,z,px,py,pz=param()
  ps=set_ps(w,x,y,z,px,py,pz)
  # Find analytical relation for precision + test X=..., erange=...
  brange=[0.0,5.0]
  bwidth=0.1
  x=np.arange(brange[0],brange[1],bwidth)
  
  b,h=ps.hist.h1('x',brange=brange,bwidth=bwidth)
  
  ih=np.interp(x,b,h) # interpolated h
  ah=sum(w) * np.exp(-x) # analytical h
  
  assert np.allclose(ih,ah,rtol=ih.std()) # TODO: ih.std() need to be changed
  
  print("max b %s"%max(b))
  print("max x %s"%max(x))
  
  assert False
  #plt.plot(x,(ih - ah)/ah)
  #input("wait")
  # test X=..., erange
  
      
def test_hist_h2():
  pass
  
def test_hist_h3():
  pass
  
def test_plot_h1():
  ps=get_ps()
  ps.plot.autoclear=False
  import matplotlib.pyplot as plt
  plt.figure(1)
  plt.clf()
  brange=[0.0,10.0]
  bwidth1=0.01
  bwidth2=0.1
  bwidth3=1.0
  x=np.arange(brange[0],brange[1],bwidth1)
  for bwidth in [bwidth1,bwidth2,bwidth3]:
    ps.plot.h1('x',bwidth=bwidth,brange=brange,label=['x','x bwidth=%.2f'%bwidth])
  plt.plot(x,sum(ps.raw.w)*np.exp(-x/1.0),label='ntot * exp(-x/1.0)')
  plt.title('h1 plot test')
  plt.legend()
  
def test_plot_h2():
  ps=get_ps()
  ps.plot.autoclear=True
  import matplotlib.pyplot as plt
  plt.figure(2)
  brange=[-10.0,10.0]
  bwidth1=1.0
  bwidth2=0.1
  ps.plot.h2('y','z',
              label=['y ($\sigma=1.0$, bwidth=%s)'%bwidth1,'z ($\sigma=2.0$, bwidth=%s)'%bwidth2,''],
              log=True,bwidth1=bwidth1,bwidth2=bwidth2,brange1=brange,brange2=brange)
  plt.title('h2 plot test')
  plt.legend()
  
def test_plot_h3():
  pass
