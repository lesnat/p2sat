#coding:utf8
import _module
import numpy as np

import unittest

class test_PhaseSpace(unittest.TestCase):
  def setUp(self):
    """
    ps=_module._PhaseSpace()
    ps.raw=ps._Raw()
    ps.hist=ps._Hist()
    ps.plot=ps._Plot()
    """
    self.ps=_module.PhaseSpaceSmilei()
    size=100000
    
    
    self.w=np.random.uniform(low=0.0,high=100.0,size=size)
    #self.w=[1.0]*size
    self.ntot=sum(self.w)
    
    self.x=np.random.exponential(scale=1.0,size=size)
    self.y=np.random.normal(loc=0.0,scale=1.0,size=size)
    self.z=np.random.normal(loc=0.0,scale=2.0,size=size)
    
    self.px=np.random.exponential(scale=1.0,size=size)
    self.py=np.random.exponential(scale=1.0,size=size)
    self.pz=np.random.exponential(scale=1.0,size=size)
    
    self.ps.raw.update(self.w,self.x,self.y,self.z,self.px,self.py,self.pz,verbose=False)
  
  def tearDown(self):
    del self.ps
    
  def test_update(self):
    print("\nLaunching test_update ...\n\n")
    np.testing.assert_almost_equal(self.w,self.ps.raw.w)
    
    np.testing.assert_almost_equal(self.x,self.ps.raw.x)
    np.testing.assert_almost_equal(self.y,self.ps.raw.y)
    np.testing.assert_almost_equal(self.z,self.ps.raw.z)
    
    np.testing.assert_almost_equal(self.px,self.ps.raw.px)
    np.testing.assert_almost_equal(self.py,self.ps.raw.py)
    np.testing.assert_almost_equal(self.pz,self.ps.raw.pz)
    
  def test_data(self):
    print("\nLaunching test_data ...\n\n")
    fname="test_export.dat"
    self.ps.raw.export_data(fname)
    #                   w   x   y   z  px  py  pz
    self.ps.raw.update([0],[0],[0],[0],[0],[0],[0])
    self.ps.raw.import_data(fname)
    
    self.test_update()
    
    
  def test_hist_hn_default(self):
    print("\nLaunching test_hist_hn_default ...\n")
    b,h=self.ps.hist.hn(['px'])
    
    b = b[0]
    
    # Test default bin length
    self.assertEqual(len(b),10)
    
    # Test default brange
    self.assertAlmostEqual(b[-1],max(self.ps.raw.x))
    self.assertAlmostEqual(b[0],min(self.ps.raw.x))
    
    # Test erange selection
    
  def test_hist_hn_erange(self):
    pass
  
  def test_hist_h1(self):
    # Find analytical relation for precision + test X=..., erange=...
    brange=[0.0,5.0]
    bwidth=0.1
    x=np.arange(brange[0],brange[1],bwidth)
    b,h=self.ps.hist.h1('x',brange=brange,bwidth=bwidth)
    ih=np.interp(x,b,h) # interpolated h
    ah=self.ntot * np.exp(-x) # analytical h
    np.testing.assert_almost_equal(ih,ah)
    #plt.plot(x,(ih - ah)/ah)
    #input("wait")
    # test X=..., erange
    
        
  def test_hist_h2(self):
    pass
    
  def test_hist_h3(self):
    pass
    
  def test_plot_h1(self):
    self.ps.plot.autoclear=False
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.clf()
    brange=[0.0,10.0]
    bwidth1=0.01
    bwidth2=0.1
    bwidth3=1.0
    x=np.arange(brange[0],brange[1],bwidth1)
    for bwidth in [bwidth1,bwidth2,bwidth3]:
      self.ps.plot.h1('x',bwidth=bwidth,brange=brange,label=['x','x bwidth=%.2f'%bwidth])
    plt.plot(x,self.ntot*np.exp(-x/1.0),label='ntot * exp(-x/1.0)')
    plt.title('h1 plot test')
    plt.legend()
    
  def test_plot_h2(self):
    self.ps.plot.autoclear=True
    import matplotlib.pyplot as plt
    plt.figure(2)
    brange=[-10.0,10.0]
    bwidth1=1.0
    bwidth2=0.1
    self.ps.plot.h2('y','z',
                    label=['y ($\sigma=1.0$, bwidth=%s)'%bwidth1,'z ($\sigma=2.0$, bwidth=%s)'%bwidth2,''],
                    log=True,bwidth1=bwidth1,bwidth2=bwidth2,brange1=brange,brange2=brange)
    plt.title('h2 plot test')
    plt.legend()
    
  def test_plot_h3(self):
    pass
    
if __name__== '__main__':
  unittest.main()
  input("Waiting ...")

