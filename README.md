# OVERVIEW

dap is a program to "Display Atoms" with "Python", based on [ASE](https://wiki.fysik.dtu.dk/ase) and [vtk](https://pypi.org/project/vtk)

# FUNCTIONALITY

## Basic Objects:
  - display atoms from any `ase.io.read` compatible format
    - atom color/radius by fixed value or based on colormap of arbitrary scalar `Atoms.arrays` property
    - atomic labels from ID# or other property (e.g. Z, species, or arbitrary `Atoms.arrays` property)
  - create bonds by cutoff or indices (optionally filtered by type)
    - bonds across PBCs displayed as two half bonds
    - multiple sets of bonds can be defined, each with own radius and color
  - display isosurfaces from volumetric data over the config unit cell (including non-orthorhombic)
    - a simple native ASCII format
    - CHGCAR using `VaspChargeDensity` from `ase.calculators.vasp `
    - WAVECAR using `Wavecar` from `pymatgen.io.vasp.outputs` (patched to read gamma-point only runs)
  - display vectors associated with each atom (like magnetic moments), either from real 3-vector field, or from scalar field (up/down, can be colored by sign)
  - display coordination polyhedra given atom type of center, and either existing bonds or cutoff to neighbors
  - display periodic images (or slice original image) along cell vectors within some (floating point) range of lattice coordinates

## Display
  - multiple configurations interpreted as trajectories (i.e. identical camera view)
  - display legend with each used atom type
  - frame label from arbitrary `Atoms.info` field
  - labels (frame and atom) do python evaluation of expressions in $(EXPR) and substitute `Atoms.info` for ${INFO} and `Atoms.arrays` for $${ARRAY}

## File Output:
  - arbitrary snapshots by rendering at integer multiples of display resolution
  - automated movie creation (built-in raw frame output, or actual movie file with ffmpeg-python package and ffmpeg executable in path)

## Other:
  - arbitrary python/ASE functions applied to atomic configurations
  - save and restore full state of display (settings, atoms, bonds, etc., as well as view orientation, but not window size) 
    when using `ase.io.write` format that supports info and arrays, e.g. extxyz.
  - command history using GNU readline 

# GUI
  - rotate, translate, zoom (zoom can also be controlled by command line)
  - pick by point or area (picked status stored as `Atoms.arrays` attribute so it is available to python functions)
  - measure picked objects (position, distance, angle), quantities displayed in _CLI_ window

****************************************************************************************************

# USAGE:
```
dap [ -g WIDTHxHEIGHT ] [ -e 'command ... [ ; command ... ] ' ] atoms_filename [ atoms_filename ... ]
```
Input atoms files can be any `ase.io.read` compatible format.

Settings commands (which are required to define atomic/bond types before the atoms files are read in) 
are read from '~/.daprc' and then `$PWD/.daprc'

Additional commands can be passed with '-e', and these include 'read commands\_filename' to read 
commands from a file.

GUI help is available by typing 'h' in the GUI window (but displayed in the CLI window).
CLI help is available with 'usage' or, for each command, with 'command -h' in the CLI winow

****************************************************************************************************

# TODO (in expected order of implementation)
##    MAYBE DONE?
  - get rid of bond\_type and surface\_type, replace with direct color/material/radius arguments to bond, volume, and polyhedra commands

##    STARTED
  - optimize rebuilding of scene to minimize cost, e.g. if only color changes are possible, don't recreate full objects unnecessarily (medium) [STARTED]

##    NOT EVEN STARTED
  - color isosurfaces by a scalar (medium)
  - measure angle for atoms, not just bonds ? (easy, but may generate too much output - add -verbose flag?)
  - display measurements in GUI ? (medium)
  - better examples (easy)
  - slices through volumetric data (maybe also voxel visualization) (depends on vtk support)
  - better integration with ASE data structures (varies)
  - read CHG as well as CHGCAR (easy)
