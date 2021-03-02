# *d*isplay *a*toms with *p*ython (dap)

Make sure the python3 and pip3 executables you want are in your path.

To install directly from this source, do (in this directory which contains the dap executable)
```
    pip3 install .
```

You can also create a wheel for for distribution with
```
    make clean
    make dist
```
and install it with
```
    pip3 install dist/dap-*-py3-none-any.whl
```
(assuming only one version's wheel file is present)

The installation process should also install prerequisites (ASE and vtk) using pip if needed.  Note that pymatgen is required for WAVECAR plotting, but that is _not_ automatically installed as a dependency - you'll need to pip install it manually.
