# -*- coding: utf-8 -*-

# https://setuptools.readthedocs.io/en/latest/setuptools.html
# https://python-packaging.readthedocs.io/en/latest/minimal.html
# https://python-guide-pt-br.readthedocs.io/en/latest/writing/structure/

from setuptools import setup

setup(name='beamFE',
      version='0.1',
      description='2D Beam finite element package for linear static and modal analysis',
      url='http://github.com/jkbgbr/beamFE2',
      author='Jakab Gabor',
      author_email='jkbgbr@gmail.com',
      license='GNU GPLv3',
      packages=['BeamFE2'],
      zip_safe=False)
