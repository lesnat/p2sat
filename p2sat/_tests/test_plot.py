#coding:utf8
import numpy as np
import matplotlib.pyplot as plt

from fixtures import *

plt.figure(figsize=(9.5,11.0))

def test_plot_h1():
  ps=get_ps()
  brange=[0.0,10.0]
  bwidth1=0.01
  bwidth2=0.1
  bwidth3=1.0
  x=np.arange(brange[0],brange[1],bwidth1)

  ps.plot.autoclear=False
  plt.subplot(311)
  for bwidth in [bwidth1,bwidth2,bwidth3]:
    ps.plot.h1('px',bwidth=bwidth,brange=brange,label=['px','bwidth=%.2f'%bwidth])
  plt.plot(x,sum(ps.raw.w)*np.exp(-x/1.0),label='ntot * exp(-x/1.0)')
  plt.title('h1,h2,h3 plot test')
  plt.legend()
  
def test_plot_h2():
  ps=get_ps()
  brange=[-10.0,10.0]
  bwidth1=1.0
  bwidth2=0.1

  ps.plot.autoclear=False
  plt.subplot(312)
  ps.plot.h2('y','z',
              label=['y ($\sigma=1.0$, bwidth=%s)'%bwidth1,'z ($\sigma=2.0$, bwidth=%s)'%bwidth2,''],
              log=True,bwidth1=bwidth1,bwidth2=bwidth2,brange1=brange,brange2=brange)
  plt.legend()
  
def test_plot_contour():
  assert False
  
def test_plot_h3():
  plt.subplot(313)
  plt.show()
  assert False


