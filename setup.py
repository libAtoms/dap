#!/usr/bin/env python3

from setuptools import setup, find_namespace_packages
import os, glob

setup(name='dap',
      version='0.9.3',
      description='Display atomic configuration',
      author='Noam Bernstein',
      author_email='noam.bernstein@nrl.navy.mil',
      packages=find_namespace_packages(where="src"),
      package_dir={"": "src"},
      package_data={'davtk.scripts': ['*'], 'davtk.examples': ['*']},
      entry_points={'console_scripts': ['dap = davtk.cli.dap:main', 'dap_ipy = davtk.cli.dap_ipy:main']},
      install_requires=['ase', 'vtk', 'PyQt6', 'numpy', 'scipy']
     )

print("")
print("See daprc for some example settings, and make your own in ~/.daprc")
print("")
print("More examples at example.* in this directory, view with")
print('     dap -e "read example.settings example.commands" example.xyz')

