# coding:utf8

import numpy as np

p2sat_path="../../"
import sys
if p2sat_path not in sys.path: sys.path.append(p2sat_path)
import p2sat

class TestData:
    def test_get_ps(self):
        ps = p2sat.ExamplePhaseSpace()
        r = ps.data.raw
        data = ps.data.get_ps()
        w,ps = data[0], data[1:]

        assert np.allclose(w,r.w)
        assert np.allclose(ps[0],r.x)
        assert np.allclose(ps[1],r.y)
        assert np.allclose(ps[2],r.z)
        assert np.allclose(ps[3],r.px)
        assert np.allclose(ps[4],r.py)
        assert np.allclose(ps[5],r.pz)
        assert np.allclose(ps[6],r.t)

    def test_get_axis(self):
        ps = p2sat.ExamplePhaseSpace()
        r = ps.data.raw

        w = ps.data.get_axis("w")
        ekin = ps.data.get_axis("ekin")

        assert np.allclose(w,r.w)
        assert np.allclose(ekin,r.ekin)

    def test_update(self):
        ps = p2sat.ExamplePhaseSpace()

        data = [np.random.rand(100) for _ in range(8)]

        ps.data.update(*data)

        for ix,ax in enumerate(ps.data.get_ps()):
            assert np.allclose(ax,data[ix])

    def test_generate(self):
        ps = p2sat.ExamplePhaseSpace()
        pass # TODO

    def test_filter_axis(self):
        ps = p2sat.ExamplePhaseSpace()
        r = ps.data.raw

        x = ps.data.filter_axis('x',select=dict(r=[None,None]))

        assert np.allclose(x,r.x)

        id = ps.data.filter_axis('id',select=dict(x=[0,None],r=[0,10],t=[100,300]))

        for i in id:
            assert r.x[i] > 0.
            assert r.r[i] > 0.
            assert r.r[i] < 10.
            assert r.t[i] > 100.
            assert r.t[i] < 300.

    def test_filter_ps(self):
        ps1 = p2sat.ExamplePhaseSpace()
        ps2 = ps1.copy()

        r1 = ps1.data.raw
        r2 = ps2.data.raw

        filter_dict = dict(x=[0,None],r=[0,10],t=[100,300])
        ps2.data.filter_ps(select=filter_dict, update=True)
        id1 = ps1.data.get_axis('id',select=filter_dict)

        assert np.allclose(r1.w[id1],  r2.w)
        assert np.allclose(r1.x[id1],  r2.x)
        assert np.allclose(r1.y[id1],  r2.y)
        assert np.allclose(r1.z[id1],  r2.z)
        assert np.allclose(r1.px[id1], r2.px)
        assert np.allclose(r1.py[id1], r2.py)
        assert np.allclose(r1.pz[id1], r2.pz)
        assert np.allclose(r1.t[id1],  r2.t)

    def test_transformate(self):
        ps1 = p2sat.ExamplePhaseSpace()
        ps2 = ps1.copy()

        r1 = ps1.data.raw
        r2 = ps2.data.raw

        ps2.data.transformate(R=(0,0,180))
        assert np.allclose(-r2.px,r1.px)
        assert np.allclose(-r2.py,r1.py)

        ps2.data.transformate(T=(-100,0,0),R=(0,0,-180),rotate_first=True)
        assert np.allclose(r2.px,r1.px)
        assert np.allclose(r2.py,r1.py)
        assert np.allclose(r2.x + 100,r1.x)

    def test_propagate(self):
        ps = p2sat.ExamplePhaseSpace()
        r = ps.data.raw

        ps.data.propagate(x=100)
        assert set(r.x) == {100}

        ps.data.propagate(t=100)
        assert set(r.t) == {100}

    def test_rescale_axis(self):
        ps1 = p2sat.ExamplePhaseSpace()
        ps2 = ps1.copy()

        ps1.data.rescale_axis("w",100.)

        assert np.allclose(ps1.data.raw.w, 100 * ps2.data.raw.w)

    def test_round_axis(self):
        ps1 = p2sat.ExamplePhaseSpace()
        ps2 = ps1.copy()

        ps1.data.round_axis("x",decimals=2)

        assert np.allclose(ps1.data.raw.x, np.around(ps2.data.raw.x,decimals=2))

    def test_deduplicate_ps(self):
        ps = p2sat.ExamplePhaseSpace()
        pass # TODO

    def test_rebin_axis(self):
        ps1 = p2sat.ExamplePhaseSpace()
        ps2 = ps1.copy()
        r1 = ps1.data.raw
        r2 = ps2.data.raw

        nbins = 10

        x = ps1.data.rebin_axis("x",nbins=nbins)

        assert len(set(x)) == nbins
        assert np.all(abs(x - r2.x) < (max(r2.x) - min(r2.x))/(nbins-1.))
        assert np.nan_to_num(sum(x)) != 0.

    def test_rebin_ps(self):
        ps = p2sat.ExamplePhaseSpace()
        pass # TODO
