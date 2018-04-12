#coding:utf8
from PhaseSpace import PhaseSpaceGeneric,PhaseSpaceSmilei,PhaseSpaceGeant4


def testing(name="all"):
  from ._tests import run
  if name=="all":
    run(name="raw")
    run(name="hist")
    run(name="plot")
    run(name="heritage")
  else:
    run(name=name)
