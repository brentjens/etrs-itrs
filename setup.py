#!/usr/bin/env python

from distutils.core import setup
import etrsitrs

setup(name='etrs-itrs',
      version=etrsitrs.__version__,
      description='A python tool to convert ETRS to ITRF coordinates',
      author='Michiel Brentjens',
      author_email='brentjens@astron.nl',
      url='',
      packages=['etrsitrs'],
      scripts=[],
      requires=['numpy', 'nose', 'coverage', 'numpydoc', 'sphinx'],
     )
