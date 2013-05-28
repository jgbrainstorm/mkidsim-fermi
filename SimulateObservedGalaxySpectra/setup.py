#!/usr/bin/env python
# encoding: utf-8
 
import sys
import os

try:
    from setuptools import setup, Extension
    setup, Extension
except ImportError:
    from distutils.core import setup, Extension
    setup, Extension

import numpy.distutils.misc_util

required = ["numpy"]

print numpy.distutils.misc_util.get_numpy_include_dirs();

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="kcorrect spokes module",
    version="0.0.1",
    author="",
    author_email="",
    url="",
    packages=["kcorrect"],
    description="python kcorrect module",
    install_requires=required,
    ext_modules=[Extension(
        "kcorrect", 
        sources = ["kcorrect/idl_k_binspec.c",  "kcorrect/k_binspec.c", "kcorrect/k_locate.c",  "kcorrect/k_midpnt.c", "kcorrect/k_polint.c",  "kcorrect/k_python.c",  "kcorrect/k_qromo.c",  "kcorrect/k_utils.c"],
	include_dirs = ["kcorrect"],
        extra_link_args = []
    )],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    classifiers=[
        #see: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: Other/Proprietary License", # TODO: this neetd to be fixed!
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Topic :: Scientific/Engineering :: Astronomy"
    ]
)
