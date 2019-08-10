from __future__ import print_function

import numpy as np, sys, re, tempfile, os

from davtk.parse_utils import ThrowingArgumentParser
try:
    from davtk.util_global import *
except ImportError:
    pass
try:
    sys.path.append('.')
    from dap_util import *
    del sys.path[-1]
except ImportError:
    pass

parsers = {}

def parse_usage(davtk_state, renderer, args):
    print("\nSETTINGS:")
    for keyword in sorted(davtk_state.settings.parsers.keys()):
        print(keyword, davtk_state.settings.parsers[keyword][1], end='')
    print("\nCOMMANDS:")
    for keyword in sorted(parsers.keys()):
        print(keyword, parsers[keyword][1], end='')
parsers["usage"] = (parse_usage, "usage: usage\n", "usage: usage\n")

def parse_help(davtk_state, renderer, args):
    for keyword in sorted(davtk_state.settings.parsers.keys()):
        print("--------------------------------------------------------------------------------")
        print(keyword, davtk_state.settings.parsers[keyword][2], end='')
    for keyword in sorted(parsers.keys()):
        print("--------------------------------------------------------------------------------")
        print(keyword, parsers[keyword][2], end='')
parsers["help"] = (parse_help, "usage: help\n", "usage: help\n")

################################################################################

parser_atom_type_field = ThrowingArgumentParser(prog="atom_type_field",description="ASE at.arrays field to use for atom type")
parser_atom_type_field.add_argument("-all_frames",action="store_true")
parser_atom_type_field.add_argument("field",type=str,help="name of field ('Z' for atomic numbers, 'species' for chemical symbols", default='Z')
def parse_atom_type_field(davtk_state, renderer, args):
    args = parser_atom_type_field.parse_args(args)
    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]

    for at in ats:
        at.info["_vtk_type_field"] = args.field
    return "cur"
parsers["atom_type_field"] = (parse_atom_type_field, parser_atom_type_field.format_usage(), parser_atom_type_field.format_help())

parser_movie = ThrowingArgumentParser(prog="movie",description="make a movie")
parser_movie.add_argument("-range",type=str,help="range of configs, in slice format start:[end+1]:[step]",default="::")
parser_movie.add_argument("-fps",type=float,help="frames per second display rate", default=10.0)
parser_movie.add_argument("-tmpdir",type=str,help="temporary directory for snapshots",default=".")
parser_movie.add_argument("-ffmpeg_args",type=str,help="other ffmpeg args",default="-b:v 10M -pix_fmt yuv420p")
parser_movie.add_argument("output_file",type=str,help="output file name")
def parse_movie(davtk_state, renderer, args):
    args = parser_movie.parse_args(args)
    m = re.search("^(\d*)(?::(\d*)(?::(\d*))?)?$", args.r)
    range_start = 0 if m.group(1) is None or len(m.group(1)) == 0 else int(m.group(1))
    range_end = len(davtk_state.at_list) if m.group(2) is None or len(m.group(2)) == 0 else int(m.group(2))
    range_interval = 1 if m.group(3) is None or len(m.group(3)) == 0 else int(m.group(3))
    frames = list(range(range_start, range_end, range_interval))
    frames.append(frames[-1])
    print(len(frames), "log",np.log10(len(frames)-1)+1)
    fmt_core = "0{}d".format(int(np.log10(len(frames)-1)+1))
    py_fmt = ".{{:{}}}.png".format(fmt_core)
    tmpfiles = []
    with tempfile.NamedTemporaryFile(dir=args.tmpdir) as fout:
        img_file_base = fout.name
        for (img_i, frame_i) in enumerate(frames):
            davtk_state.show_frame(frame_i = frame_i)
            tmpfiles.append(img_file_base+py_fmt.format(img_i))
            davtk_state.snapshot(tmpfiles[-1], 1)
    print("ffmpeg -i {}.%{}.png -r {} -b:v 10M -pix_fmt yuv420p {}".format(img_file_base, fmt_core, args.fps, args.output_file))
    os.system("ffmpeg -i {}.%{}.png -r {} {} {}".format(img_file_base, fmt_core, args.fps, args.ffmpeg_args, args.output_file))
    for f in tmpfiles:
        os.remove(f)
    return None
parsers["movie"] = (parse_movie, parser_movie.format_usage(), parser_movie.format_help())

parser_go = ThrowingArgumentParser(prog="go",description="go to a particular frame")
parser_go.add_argument("n",type=int,help="number of frames to change")
def parse_go(davtk_state, renderer, args):
    args = parser_go.parse_args(args)
    davtk_state.show_frame(frame_i = args.n)
    return None
parsers["go"] = (parse_go, parser_go.format_usage(), parser_go.format_help())

parser_next = ThrowingArgumentParser(prog="next",description="go forward a number of frames")
parser_next.add_argument("n",type=int,nargs='?',default=0,help="number of frames to change (default set by 'step' command)")
def parse_next(davtk_state, renderer, args):
    args = parser_next.parse_args(args)
    davtk_state.show_frame(dframe = args.n if args.n > 0 else davtk_state.settings["frame_step"])
    return None
parsers["next"] = (parse_next, parser_next.format_usage(), parser_next.format_help())

parser_prev = ThrowingArgumentParser(prog="prev", description="go back a number of frames")
parser_prev.add_argument("n",type=int,nargs='?',default=0, help="number of frames to change (default set by 'step' command)")
def parse_prev(davtk_state, renderer, args):
    args = parser_prev.parse_args(args)
    davtk_state.show_frame(dframe = -args.n if args.n > 0 else -davtk_state.settings["frame_step"])
    return None
parsers["prev"] = (parse_prev, parser_prev.format_usage(), parser_prev.format_help())

parser_unpick = ThrowingArgumentParser(prog="unpick", description="unpick all picked atoms")
parser_unpick.add_argument("-all_frames",action="store_true")
parser_unpick.add_argument("-only",nargs='+',choices=["atoms","bonds"],action="store", default=None)
def parse_unpick(davtk_state, renderer, args):
    args = parser_unpick.parse_args(args)
    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]
    for at in ats:
        if args.only is None or "atoms" in args.only:
            at.arrays["_vtk_picked"][:] = False
        if args.only is None or "bonds" in args.only:
            if hasattr(at, "bonds"):
                for i_at in range(len(at)):
                    for b in at.bonds[i_at]:
                        b["picked"] = False
    return "cur"
parsers["unpick"] = (parse_unpick, parser_unpick.format_usage(), parser_unpick.format_help())

parser_pick = ThrowingArgumentParser(prog="pick", description="pick atom(s) by ID")
parser_pick.add_argument("-all_frames",action="store_true")
parser_pick.add_argument("n",type=int,nargs='+')
def parse_pick(davtk_state, renderer, args):
    args = parser_pick.parse_args(args)
    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]
    for at in ats:
        at.arrays["_vtk_picked"][args.n] = True
    return "all"
parsers["pick"] = (parse_pick, parser_pick.format_usage(), parser_pick.format_help())

parser_delete = ThrowingArgumentParser(prog="delete",description="delete objects (picked by default)")
parser_delete.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_delete.add_argument("-atoms",type=int,nargs='+',help="delete by these indices ")
parser_delete.add_argument("-bonds",action="store_true",help="delete all bonds")
def parse_delete(davtk_state, renderer, args):
    args = parser_delete.parse_args(args)
    if args.all_frames:
        frame_list=None
    else:
        frame_list="cur"
    if args.atoms is None:
        davtk_state.delete(atoms="picked", bonds="picked", frames=frame_list)
    else:
        davtk_state.delete(atoms=args.atoms, frames=frame_list)

    if args.bonds:
        davtk_state.delete(bonds="all", frames=frame_list)
    return None
parsers["delete"] = (parse_delete, parser_delete.format_usage(), parser_delete.format_help())

parser_dup = ThrowingArgumentParser(prog="dup",description="duplicate cell")
parser_dup.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_dup.add_argument("n",type=int,nargs='+', help="number of images to set (scalar or 3-vector)")
def parse_dup(davtk_state, renderer, args):
    args = parser_dup.parse_args(args)
    if args.all_frames:
        frame_list=None
    else:
        frame_list="cur"
    if len(args.n) == 1:
        args.n *= 3
    elif len(args.n) != 3:
        raise ValueError("wrong number of elements (not 1 or 3) in n "+str(args.n))
    davtk_state.duplicate(args.n, frames=frame_list)
    return None
parsers["dup"] = (parse_dup, parser_dup.format_usage(), parser_dup.format_help())

parser_images = ThrowingArgumentParser(prog="images",description="show images of cell")
parser_images.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_images.add_argument("n",type=float,nargs='+', help="number of images to set (scalar or 3-vector)")
def parse_images(davtk_state, renderer, args):
    args = parser_images.parse_args(args)

    if len(args.n) == 1:
        args.n *= 3
    elif len(args.n) != 3:
        raise ValueError("'"+str(args.n)+" not scalar or 3-vector")

    if args.all_frames:
        for at in davtk_state.at_list:
            at.info["_vtk_images"] = args.n
    else:
        davtk_state.cur_at().info["_vtk_images"] = args.n

    davtk_state.show_frame(dframe=0)

    return None
parsers["images"] = (parse_images, parser_images.format_usage(), parser_images.format_help())

parser_bond = ThrowingArgumentParser(prog="bond",description="Create bonds")
parser_bond.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_bond.add_argument("-name",type=str,help="name of bond type", default=None)
parser_bond.add_argument("-T",type=str,help="Restrict one member to given atom_type", default="*")
parser_bond.add_argument("-T2",type=str,help="Restrict other member to given atom_type", default="*")
grp = parser_bond.add_mutually_exclusive_group()
grp.add_argument("-picked", action='store_true', help="bond a pair of picked atoms with MIC")
grp.add_argument("-n", action='store', type=int, nargs=2, help="bond a pair of atoms with given IDs", default=None)
grp.add_argument("-rcut", dest="cutoff", action='store', type=float, nargs='+', help="bond atoms with max (one argument) or min-max (two arguments) cutoffs", default=None)
def parse_bond(davtk_state, renderer, args):
    args = parser_bond.parse_args(args)

    if args.cutoff is not None and len(args.cutoff) > 2:
        raise ValueError("bond got -rcut with {} > 2 values".format(len(args.cutoff)))

    if args.all_frames:
        frames = None
    else:
        frames = "cur"

    if args.picked:
        davtk_state.bond(args.name, args.T, args.T2, "picked", frames)
    elif args.n:
        davtk_state.bond(args.name, args.T, args.T2, ("n", args.n), frames)
    elif args.cutoff:
        davtk_state.bond(args.name, args.T, args.T2, ("cutoff", args.cutoff), frames)
    else:
        davtk_state.bond(args.name, args.T, args.T2, ("cutoff", None), frames)
    return None
parsers["bond"] = (parse_bond, parser_bond.format_usage(), parser_bond.format_help())

parser_snapshot = ThrowingArgumentParser(prog="snapshot",description="write snapshot")
parser_snapshot.add_argument("-mag",type=int,help="magnification", default=1)
parser_snapshot.add_argument("file",type=str,help="filename")
def parse_snapshot(davtk_state, renderer, args):
    args = parser_snapshot.parse_args(args)
    davtk_state.snapshot(args.file, args.mag)
parsers["snapshot"] = (parse_snapshot, parser_snapshot.format_usage(), parser_snapshot.format_help())

parser_X = ThrowingArgumentParser(prog="X",description="execute python code (Atoms object available as 'atoms', DavTKSettings as 'settings')")
parser_X.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_X.add_argument("-global_context","-g",action="store_true",help="run in global scope (for 'import', e.g.)")
parser_X.add_argument("args",type=str,nargs='+',help="X command line words")
def parse_X(davtk_state, renderer, args):
    args = parser_X.parse_args(args)
    ase_command = " ".join(args.args)
    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]
    settings = davtk_state.settings
    for atoms in ats:
        globals()["atoms"] = atoms
        if args.global_context:
            exec(ase_command) in globals(), globals()
        else:
            exec(ase_command)
    return "all"
parsers["X"] = (parse_X, parser_X.format_usage(), parser_X.format_help())

parser_read = ThrowingArgumentParser(prog="read",description="read commands from file(s)")
parser_read.add_argument("filename",type=str,nargs='+',help="filenames")
def parse_read(davtk_state, renderer, args):
    args = parser_read.parse_args(args)
    for f in args.filename:
        with open(f) as fin:
            for l in fin.readlines():
                parse_line(l, davtk_state.settings, davtk_state, davtk_state.renderer)
parsers["read"] = (parse_read, parser_read.format_usage(), parser_read.format_help())

parser_label = ThrowingArgumentParser(prog="label",description="enable/disable labels (toggle by default)")
parser_label.add_argument("-all_frames",action="store_true",help="apply to all frames")
group = parser_label.add_mutually_exclusive_group()
group.add_argument("-on",action='store_true',help="turn on")
group.add_argument("-off",action='store_true',help="turn off")
def parse_label(davtk_state, renderer, args):
    args = parser_label.parse_args(args)
    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]

    for at in ats:
        if args.on:
            at.info["_vtk_show_labels"] = True
        elif args.off:
            at.info["_vtk_show_labels"] = False
        else:
            at.info["_vtk_show_labels"] = not at.info["_vtk_show_labels"] 
    return refresh
parsers["label"] = (parse_label, parser_label.format_usage(), parser_label.format_help())

parser_measure = ThrowingArgumentParser(prog="measure",description="measure some quantities for picked objects (default) or listed atoms")
parser_measure.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_measure.add_argument("-n",action="store",type=int,nargs='+',help="ID numbers of atoms to measure", default = None)
def parse_measure(davtk_state, renderer, args):
    args = parser_measure.parse_args(args)
    if args.all_frames:
        frames = range(len(davtk_state.at_list))
    else:
        frames = [davtk_state.cur_frame]

    for frame_i in frames:
        print("Frame:",frame_i)
        davtk_state.measure(args.n, frame_i)
    return None
parsers["measure"] = (parse_measure, parser_measure.format_usage(), parser_measure.format_help())

################################################################################

def parse_line(line, settings, state, renderer=None):
    if len(line.strip()) > 0:
        args = line.split()

        matches_settings = [ k for k in settings.parsers.keys() if k.startswith(args[0]) ]
        matches_cmds = [ k for k in parsers.keys() if k.startswith(args[0]) ]

        if len(matches_settings+matches_cmds) == 0: # no match
            raise ValueError("Unknown command '{}'\n".format(args[0]))

        if len(matches_settings+matches_cmds) == 1: # unique match
            if len(matches_settings) == 1:
                return settings.parsers[matches_settings[0]][0](args[1:])
            else: # must be matches_cmds
                return parsers[matches_cmds[0]][0](state, renderer, args[1:])

        if args[0] in matches_settings:
            return settings.parsers[args[0]][0](args[1:])
        if args[0] in matches_cmds:
            return parsers[args[0]][0](state, renderer, args[1:])

        raise ValueError("Ambiguous command '{}', matches {}".format(args[0], matches_settings+matches_cmds))
    else:
        return None

def parse_file(filename, settings, state=None):
    with open(filename) as fin:
        for l in fin.readlines():
            if l.startswith("#"):
                continue
            refresh = parse_line(l, settings, state)
            if state is not None:
                state.update(refresh)
