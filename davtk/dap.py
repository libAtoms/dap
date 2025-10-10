import argparse
from .__version__ import __version__

cli_parse = argparse.ArgumentParser()
cli_parse.add_argument("--geometry","-g",type=str,default="1000x1000")
cli_parse.add_argument("--execute_commands","-e",action='append',type=str, help="commands to execute (';' separated)", default=[])
cli_parse.add_argument("--format","-f",type=str, nargs='+', help="formats for each atoms file", default = None)
cli_parse.add_argument("--version","-v", action='version', version='%(prog)s ' + __version__)
cli_parse.add_argument("atoms_files",nargs="+",type=str,help="atoms files in any ase.io.read format, optionally including '@START:END:STEP' interval")
args = cli_parse.parse_args()

import sys, re
import ase, ase.io
from davtk.viewer import Viewer

# read atoms from atoms_files
if args.format is None:
    args.format = [ None ] * len(args.atoms_files)
else:
    if len(args.format) != len(args.atoms_files):
        raise ValueError("Got {} formats, but different number of atoms files {}".format(len(args.format), len(args.atoms_files)))
at_list = []
for (filename, fmt) in zip(args.atoms_files, args.format):
    if "@" in filename:
        r = None
    else:
        r = ":"
    ats = ase.io.read(filename, r, format=fmt)
    if isinstance(ats, ase.Atoms):
        ats = [ ats ]
    for at in ats:
        if not "_vtk_filename" in at.info:
            at.info["_vtk_filename"] = filename
        at_list.append(at)

m = re.search("^(\d+)x(\d+)$", args.geometry)
if not m:
    raise ValueError("-geometry '{}' not NxM".format(args.geometry))
win_size = (int(m.group(1)), int(m.group(2)))
viewer = Viewer(at_list, win_size, "dap " + " ".join(args.atoms_files), args.execute_commands)
