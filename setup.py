#!/usr/bin/env python

from distutils.core import setup

setup(name='dap',
      version='0.8',
      description='Display atomic configuration',
      author='Noam Bernstein',
      author_email='noam.bernstein@nrl.navy.mil',
      py_modules=[ 'davtk_interactors',
                   'davtk_parse',
                   'davtk_parse_interactive',
                   'davtk_parse_utils',
                   'davtk_settings',
                   'davtk_state',
                   'davtk_util_global'],
     data_files=[ ('config', ['davtkrc']) ],
     scripts=['dap']
     )
