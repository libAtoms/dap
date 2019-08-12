#!/usr/bin/env python

from setuptools import setup, find_packages
import os

setup(name='dap',
      version='0.8.1',
      description='Display atomic configuration',
      author='Noam Bernstein',
      author_email='noam.bernstein@nrl.navy.mil',
      packages=find_packages(),
      data_files=[('config',['daprc', 'README'])],
      scripts=['dap'],
      install_requires=['ase','vtk']
     )

print("")
print("Don't forget to check out README")
print("See daprc (in this directory) for example settings, and/or copy it:")
print("     cp "+os.path.join(os.environ["PWD"],"daprc")+" ~/.daprc")
