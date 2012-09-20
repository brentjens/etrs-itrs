#!/usr/bin/env python

from distutils.core import setup

setup(name='etrsitrs',
      version='0.1',
      description='A python tool to convert ETRS to ITRF coordinates',
      author='Michiel Brentjens',
      author_email='brentjens@astron.nl',
      url='',
      packages=['etrsitrs'],
      scripts=[],
      requires=['numpy', 'nose', 'coverage', 'IPython(>=0.11)', 'numpydoc', 'sphinx'],
     )
