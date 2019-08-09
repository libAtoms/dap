#!/usr/bin/env python

from setuptools import setup, find_packages
import os

setup(name='dap',
      version='0.8',
      description='Display atomic configuration',
      author='Noam Bernstein',
      author_email='noam.bernstein@nrl.navy.mil',
      packages=find_packages(),
      data_files=[('config',['daprc'])],
      scripts=['dap'],
      install_requires=['ase','vtk']
     )

print("")
print("Don't forget to")
print("     cp daprc ~/.daprc")
