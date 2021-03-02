#!/usr/bin/env python3

from setuptools import setup, find_packages
import os, glob

setup(name='dap',
      version='0.8.20',
      description='Display atomic configuration',
      author='Noam Bernstein',
      author_email='noam.bernstein@nrl.navy.mil',
      packages=find_packages(),
      data_files=[('config',['daprc'] + glob.glob('example.*'))],
      scripts=['dap'],
      install_requires=['ase','vtk']
     )

print("")
print("See daprc for some example settings, and make your own in ~/.daprc")
print("")
print("More examples at example.* in this directory, view with")
print('     dap -e "read example.settings example.commands" example.xyz')

