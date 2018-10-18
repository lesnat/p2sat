#!/usr/bin/env python
#coding:utf8

from setuptools import setup

import p2sat

setup(
    name='p2sat',
    version=p2sat.__version__,
    packages=['p2sat'],
    author="lesnat",
    description="Particle Phase Space Analysis Toolkit",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    # install_requires= ['numpy','matplotlib'],
    include_package_data=True,
    url='https://github.com/lesnat/p2sat',
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "Natural Language :: English",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics"
    ]
)
