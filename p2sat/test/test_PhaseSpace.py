# coding:utf8

import numpy as np

p2sat_path="../../"
import sys
if p2sat_path not in sys.path: sys.path.append(p2sat_path)
import p2sat

class TestPhaseSpace:
    def test_import(self):
        p2sat.ExamplePhaseSpace()

    def test_particle(self):
        p2sat.PhaseSpace(particle="e-")
        p2sat.PhaseSpace(particle="electron")
        p2sat.PhaseSpace(particle="e+")
        p2sat.PhaseSpace(particle="positron")
        p2sat.PhaseSpace(particle="g")
        p2sat.PhaseSpace(particle="gamma")
        p2sat.PhaseSpace(particle="mu-")
        p2sat.PhaseSpace(particle="muon-")
        p2sat.PhaseSpace(particle="mu+")
        p2sat.PhaseSpace(particle="muon+")
        p2sat.PhaseSpace(particle="p")
        p2sat.PhaseSpace(particle="proton")
        p2sat.PhaseSpace(particle="n")
        p2sat.PhaseSpace(particle="neutron")

    def test_add(self):
        ps1 = p2sat.ExamplePhaseSpace()
        ps2 = ps1 + ps1

        assert len(ps2) == 2 * len(ps1)

        i = len(ps1)

        assert np.allclose(ps1.data.raw.w,  ps2.data.raw.w[:i])
        assert np.allclose(ps1.data.raw.x,  ps2.data.raw.x[:i])
        assert np.allclose(ps1.data.raw.y,  ps2.data.raw.y[:i])
        assert np.allclose(ps1.data.raw.z,  ps2.data.raw.z[:i])
        assert np.allclose(ps1.data.raw.px, ps2.data.raw.px[:i])
        assert np.allclose(ps1.data.raw.py, ps2.data.raw.py[:i])
        assert np.allclose(ps1.data.raw.pz, ps2.data.raw.pz[:i])
        assert np.allclose(ps1.data.raw.t,  ps2.data.raw.t[:i])

        assert np.allclose(ps1.data.raw.w,  ps2.data.raw.w[i:])
        assert np.allclose(ps1.data.raw.x,  ps2.data.raw.x[i:])
        assert np.allclose(ps1.data.raw.y,  ps2.data.raw.y[i:])
        assert np.allclose(ps1.data.raw.z,  ps2.data.raw.z[i:])
        assert np.allclose(ps1.data.raw.px, ps2.data.raw.px[i:])
        assert np.allclose(ps1.data.raw.py, ps2.data.raw.py[i:])
        assert np.allclose(ps1.data.raw.pz, ps2.data.raw.pz[i:])
        assert np.allclose(ps1.data.raw.t,  ps2.data.raw.t[i:])


    def test_copy(self):
        ps1 = p2sat.ExamplePhaseSpace()
        ps2 = ps1.copy()

        assert ps1 is not ps2

        assert np.allclose(ps1.data.raw.w,  ps2.data.raw.w)
        assert np.allclose(ps1.data.raw.x,  ps2.data.raw.x)
        assert np.allclose(ps1.data.raw.y,  ps2.data.raw.y)
        assert np.allclose(ps1.data.raw.z,  ps2.data.raw.z)
        assert np.allclose(ps1.data.raw.px, ps2.data.raw.px)
        assert np.allclose(ps1.data.raw.py, ps2.data.raw.py)
        assert np.allclose(ps1.data.raw.pz, ps2.data.raw.pz)
        assert np.allclose(ps1.data.raw.t,  ps2.data.raw.t)

    def test_len(self):
        ps = p2sat.ExamplePhaseSpace()
        assert len(ps) == len(ps.data.raw.w)
