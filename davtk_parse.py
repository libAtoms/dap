import numpy as np, sys
from davtk_parse_utils import ThrowingArgumentParser
from davtk_settings import UnknownSettingsKeywordError
try:
    from davtk_util_global import *
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
    for keyword in sorted(parsers.keys()):
        print keyword, parsers[keyword][1],
parsers["usage"] = (parse_usage, "usage: usage\n", "usage: usage\n")

def parse_help(davtk_state, renderer, args):
    for keyword in sorted(parsers.keys()):
        print "--------------------------------------------------------------------------------"
        print keyword, parsers[keyword][2],
parsers["help"] = (parse_help, "usage: help\n", "usage: help\n")

################################################################################

parser_next = ThrowingArgumentParser(prog="next",description="go forward a number of frames")
parser_next.add_argument("-n",type=int,default=1,help="number of frames to change")
def parse_next(davtk_state, renderer, args):
    args = parser_next.parse_args(args)
    davtk_state.show_frame(dframe = args.n, renderer=renderer)
    return None
parsers["next"] = (parse_next, parser_next.format_usage(), parser_next.format_help())

parser_prev = ThrowingArgumentParser(prog="prev", description="go back a number of frames")
parser_prev.add_argument("-n",type=int,default=1, help="number of frames to change")
def parse_prev(davtk_state, renderer, args):
    args = parser_prev.parse_args(args)
    davtk_state.show_frame(dframe = -args.n, renderer=renderer)
    return None
parsers["prev"] = (parse_prev, parser_prev.format_usage(), parser_prev.format_help())

parser_unpick = ThrowingArgumentParser(prog="unpick", description="unpick all picked atoms")
parser_unpick.add_argument("-all_frames",action="store_true")
def parse_unpick(davtk_state, renderer, args):
    args = parser_unpick.parse_args(args)
    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]
    for at in ats:
        at.arrays["_vtk_picked"][:] = False
        if hasattr(at, "bonds"):
            for i_at in range(len(at)):
                for b in at.bonds[i_at]:
                    b["picked"] = False
    return "all"
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

parser_bond = ThrowingArgumentParser(prog="bond",description="create bonds")
parser_bond.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_bond.add_argument("-name",type=str,help="name of bond type", default=None)
parser_bond.add_argument("a",type=float,nargs='*', help="scalar float cutoff (overriding atom type) or pair of atom index ints")
def parse_bond(davtk_state, renderer, args):
    args = parser_bond.parse_args(args)

    if args.all_frames:
        frames = None
    else:
        frames = "cur"

    if len(args.a) == 0:
        davtk_state.bond("auto_cutoff", args.name, frames)
    elif len(args.a) == 1:
        davtk_state.bond(args.a[0], args.name, frames)
    elif len(args.a) == 2:
        davtk_state.bond((int(args.a[0]),int(args.a[1])), args.name, frames)
    else:
        raise ValueError("bond got more than 2 values (1 is cutoff, 2 is index pair)")

    davtk_state.show_frame(dframe=0)
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
        refresh = "all"
    else:
        ats = [davtk_state.cur_at()]
        refresh = "cur"

    for at in ats:
        if args.on:
            at.info["_vtk_show_labels"] = True
        elif args.off:
            at.info["_vtk_show_labels"] = False
        else:
            at.info["_vtk_show_labels"] = not at.info["_vtk_show_labels"] 
    return "cur"

parsers["label"] = (parse_label, parser_label.format_usage(), parser_label.format_help())

################################################################################

def parse_line(line, settings, state, renderer=None):
    if len(line.strip()) > 0:
        try:
            return settings.parse_line(line)
        except UnknownSettingsKeywordError:
            pass

        args = line.split()
        return parsers[args[0]][0](state, renderer, args[1:])
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
