#coding:utf8

def run(name):
  import pytest
  pytest.main('/mnt/local/esnault/Modules/p2sat/_tests/test_%s.py'%name)
