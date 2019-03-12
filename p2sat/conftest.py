# coding:utf8

"""
File to make doctests with pytest.
"""

p2sat_path="../"
import sys
if p2sat_path not in sys.path: sys.path.append(p2sat_path)
import p2sat

import pytest

@pytest.fixture(autouse=True)
def add_p2sat(doctest_namespace):
    doctest_namespace['ExamplePhaseSpace'] = p2sat.ExamplePhaseSpace
