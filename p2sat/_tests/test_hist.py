#coding:utf8
import numpy as np

from fixtures import *

def test_hn_default_bins():
  ps=get_ps()
  b,h=ps.hist.hn(['px'])
  
  print('b = %s'%b)
  print('h = %s'%h)

  b = b[0]
  
  # Test default bin length
  assert len(b) == 11 # ???
  
  # Test default brange
  assert np.isclose(b[0],min(ps.raw.px))
  assert np.isclose(b[-1],max(ps.raw.px)) # ???
  
def test_hn_change_bins():
  w,x,y,z,px,py,pz=param(size=5)
  w = np.array([1.0]*5)
  px= np.array([7,4,8,8,4])
  py= np.array([3,4,6,1,4])
  ps=set_ps(w,x,y,z,px,py,pz)
  
  for bwidth in [1.0,0.5,0.2]:
    for brange1 in [[0.0,8.0],[1.0,9.0],[0.0,10.0]]:
      brange2=[0.0,8.0]
      b,h=ps.hist.hn(['px','py'],bwidth=[bwidth,bwidth],brange=[brange1,brange2])
      b1=b[0][1:] ; b2=b[1][1:] # ???
      print("bwidth = %s"%bwidth)
      print("brange1 = %s"%brange1)
      print("brange2 = %s"%brange2)
      print("b1 = %s"%b1)
      print("b2 = %s"%b2)
      print("h = \n%s"%h)
      
      assert h[b1==7.,b2==3.]==[1.0]
      assert h[b1==4.,b2==4.]==[2.0]
      assert h[b1==8.,b2==6.]==[1.0]
      assert h[b1==8.,b2==1.]==[1.0]
      
      assert h[b1==1.,b2==1.]==[0.0]
      assert h[b1==7.,b2==9.]==[0.0]
      assert h[b1==4.,b2==7.]==[0.0]

def test_hn_select():
  w,x,y,z,px,py,pz=param(size=10)
  w = np.array([1.0]*10)
  x = np.array([5,5,5,8,7,7,1,3,8,5])
  px= np.array([7,3,4,5,8,7,8,4,8,4])
  py= np.array([3,8,7,6,4,6,7,8,1,7])
  ps=set_ps(w,x,y,z,px,py,pz)
  
  bwidth=1.0
  brange=[1.0,9.0]
  for X in [5,5.0,[4.9,5.1]]:
    b,h=ps.hist.hn(['px','py'],bwidth=[bwidth,bwidth],brange=[brange,brange],select={'x':X})
    b1=b[0][:-1] ; b2=b[1][:-1] # ???
    print("b1 = %s"%b1)
    print("b2 = %s"%b2)
    print("h = \n%s"%h)
    
    assert h[b1==7.,b2==3.]==[1.0]
    assert h[b1==3.,b2==8.]==[1.0]
    assert h[b1==4.,b2==7.]==[2.0]
    
    assert h[b1==8.,b2==4.]==[0.0]
    assert h[b1==4.,b2==8.]==[0.0]
    assert h[b1==8.,b2==1.]==[0.0]
    


def test_hn_wnorm():
  ps=get_ps()
  bwidth=0.5
  assert False

def test_h1():
  w,x,y,z,px,py,pz=param()
  ps=set_ps(w,x,y,z,px,py,pz)
  # Find analytical relation for precision + test X=..., erange=...
  brange=[0.0,20.0]
  bwidth=0.5
  x=np.arange(brange[0],brange[1],bwidth)
  
  b,h=ps.hist.h1('px',brange=brange,bwidth=bwidth)
  
  ih=np.interp(x,b,h) # interpolated h
  ah=sum(w) * np.exp(-x) # analytical h
  
  assert np.allclose(ih,ah,rtol=ih.std()) # TODO: ih.std() need to be changed
  
  
  #plt.plot(x,(ih - ah)/ah)
  #input("wait")
  # test X=..., erange
  assert False
      
def test_h2():
  assert False
  
def test_h3():
  assert False


