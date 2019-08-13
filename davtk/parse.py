from __future__ import print_function

import numpy as np, sys, re, tempfile, os
import ase.io
from ase.calculators.vasp import VaspChargeDensity
import ffmpeg

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
parser_write_state = ThrowingArgumentParser(prog="write_state",description="write state of dap, including atomic configs, settings, and view")
parser_write_state.add_argument("-cur_frame_only",action="store_true")
parser_write_state.add_argument("filename",action="store",type=str,help="name of file to save into")
def parse_write_state(davtk_state, renderer, args):
    args = parser_write_state.parse_args(args)

    with open(args.filename+".settings","w") as fout:
        davtk_state.settings.write(fout)
        fout.write("restore_view -data "+" ".join([str(v) for v in davtk_state.get_view()]))

    with open(args.filename,"w") as fout:
        if args.cur_frame_only:
            ats = [davtk_state.cur_at()]
        else:
            ats = davtk_state.at_list

        davtk_state.prep_for_atoms_write(ats)

        ase.io.write(fout, ats, format=os.path.splitext(args.filename)[1].replace(".",""))
    print("Settings written to '{0}.settings', read back with\n    dap -e \"read {0}.settings\" {0}".format(args.filename))
    return None
parsers["write_state"] = (parse_write_state, parser_write_state.format_usage(), parser_write_state.format_help())

parser_restore_view = ThrowingArgumentParser(prog="restore_view",description="use stored viewing transform in ASE atoms object")
group = parser_restore_view.add_mutually_exclusive_group()
group.add_argument("-name",action="store",type=str,help="name of info field to restore view from", default="_vtk_view")
group.add_argument("-data",type=float,nargs='+',action="store",help="list of values defining view")
parser_restore_view.add_argument("-in_config",action="store_true",help="get data from configuration info dict")
def parse_restore_view(davtk_state, renderer, args):
    args = parser_restore_view.parse_args(args)

    if args.data:
        davtk_state.restore_view(args.data)
    else: # by name
        if args.in_config:
            at = davtk_state.cur_at()
            if args.name in at.info:
                davtk_state.restore_view(at.info["_vtk_view_"+args.name])
            else:
                raise ValueError("restore_view info field '{}' not found".format(args.name))
        else:
            davtk_state.restore_view(davtk_state.saved_views[args.name])
    return None
parsers["restore_view"] = (parse_restore_view, parser_restore_view.format_usage(), parser_restore_view.format_help(), None)

parser_save_view = ThrowingArgumentParser(prog="save_view",description="store viewing transform in ASE atoms object")
parser_save_view.add_argument("-all_frames",action="store_true")
parser_save_view.add_argument("-in_config",action="store_true",help="save data in configuration info dict")
parser_save_view.add_argument("-name",action="store",type=str,help="name of info field to save view into", default="_vtk_view")
def parse_save_view(davtk_state, renderer, args):
    args = parser_save_view.parse_args(args)
    view = davtk_state.get_view()
    if args.in_config:
        if args.all_frames:
            ats = davtk_state.at_list
        else:
            ats = [davtk_state.cur_at()]

        for at in ats:
            at.info["_vtk_view_"+args.name] = view
    else:
        davtk_state.saved_views[args.name] = view

    return None
parsers["save_view"] = (parse_save_view, parser_save_view.format_usage(), parser_save_view.format_help(), None)

parser_movie = ThrowingArgumentParser(prog="movie",description="make a movie")
parser_movie.add_argument("-range",type=str,help="range of configs, in slice format start:[end+1]:[step]",default="::")
parser_movie.add_argument("-framerate",type=float,help="frames per second display rate", default=10.0)
parser_movie.add_argument("-tmpdir",type=str,help="temporary directory for snapshots",default=".")
parser_movie.add_argument("-ffmpeg_args",type=str,help="other ffmpeg args",default="-b:v 10M -pix_fmt yuv420p")
parser_movie.add_argument("output_file",type=str,help="output file name")
def parse_movie(davtk_state, renderer, args):
    args = parser_movie.parse_args(args)
    m = re.search("^(\d*)(?::(\d*)(?::(\d*))?)?$", args.range)
    range_start = 0 if m.group(1) is None or len(m.group(1)) == 0 else int(m.group(1))
    range_end = len(davtk_state.at_list) if m.group(2) is None or len(m.group(2)) == 0 else int(m.group(2))
    range_interval = 1 if m.group(3) is None or len(m.group(3)) == 0 else int(m.group(3))
    frames = list(range(range_start, range_end, range_interval))

    for frame_i in frames:
        davtk_state.show_frame(frame_i)
        data = davtk_state.snapshot().astype(np.uint8)
        data = np.ascontiguousarray(np.flip(data, axis=0))

        if frame_i == frames[0]:
            # start ffmpeg
            process = (
                ffmpeg
                .input('pipe:', format='rawvideo', pix_fmt='rgb24', s='{}x{}'.format(data.shape[1],data.shape[0]))
                .output(args.output_file, pix_fmt='yuv420p', framerate=args.framerate)
                .overwrite_output()
                .run_async(pipe_stdin=True)
            )

        process.stdin.write(data)
    process.stdin.close()
    process.wait()

    # frames.append(frames[-1])
    # fmt_core = "0{}d".format(int(np.log10(len(frames)-1)+1))
    # py_fmt = ".{{:{}}}.png".format(fmt_core)
    # tmpfiles = []
    # with tempfile.NamedTemporaryFile(dir=args.tmpdir) as fout:
        # img_file_base = fout.name
        # for (img_i, frame_i) in enumerate(frames):
            # davtk_state.show_frame(frame_i = frame_i)
            # tmpfiles.append(img_file_base+py_fmt.format(img_i))
            # davtk_state.snapshot(tmpfiles[-1], 1)
    # print    ("ffmpeg -i {}.%{}.png -r {} {} {}".format(img_file_base, fmt_core, args.framerate, args.ffmpeg_args, "alt_"+args.output_file))
    # os.system("ffmpeg -i {}.%{}.png -r {} {} {}".format(img_file_base, fmt_core, args.framerate, args.ffmpeg_args, "alt_"+args.output_file))
    # for f in tmpfiles:
        # os.remove(f)

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

    got_arg = False
    if args.atoms is not None:
        davtk_state.delete(atoms=args.atoms, frames=frame_list)
        got_arg = True
    if args.bonds:
        davtk_state.delete(bonds="all", frames=frame_list)
        got_arg = True
    if not got_arg:
        davtk_state.delete(atoms="picked", bonds="picked", frames=frame_list)

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
    import ase
    from ase.io import write
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
        parse_file(f, davtk_state.settings, davtk_state)
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
    return "cur"
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

parser_volume = ThrowingArgumentParser(prog="volume",description="read volumetric data from file")
parser_volume.add_argument("filename",nargs='?',help="file to read from")
parser_volume.add_argument("-name",type=str, help="name to assign", default=None)
group = parser_volume.add_mutually_exclusive_group()
group.add_argument("-delete", action='store', metavar="NAME", help="name of volume to delete", default=None)
group.add_argument("-list", action='store_true', help="list existing volume representations")
group.add_argument("-isosurface",type=float,nargs=5,action='append',metavar=("THRESHOLD","R","G","B","OPACITY"), help="isosurface threshold, color, and opacity")
group.add_argument("-volumetric",type=float,nargs=4,action='append',metavar=("SCALE","R","G","B"), help="volumetric value_to_opacity_factor and color")
def parse_volume(davtk_state, renderer, args):
    args = parser_volume.parse_args(args)
    if args.delete is not None:
        davtk_state.delete_volume_rep(args.delete)
    elif args.list:
        print("Defined volume names: "+" ".join(davtk_state.cur_at().volume_reps.keys()))
    else:
        if args.filename is None:
            raise ValueError("volume got no filename")
        if args.name is None:
            args.name = args.filename

        if args.filename.endswith("CHGCAR"):
            chgcar = VaspChargeDensity(args.filename)
            data = chgcar.chg[0]
        else:
            with open(args.filename) as fin:
                extents = [int(i) for i in fin.readline().rstrip().split()]
                if len(extents) != 3:
                    raise ValueError("Got bad number of extents {} != 3 on first line of '{}'".format(len(full_extents), filename))

                # order of indices in data is being reversed
                data = np.zeros(extents[::-1])
                for l in fin:
                    (i0, i1, i2, v) = l.rstrip().split()
                    data[int(i2),int(i1),int(i0)] = float(v)

        if args.isosurface:
            for params in args.isosurface:
                davtk_state.add_volume_rep(args.name, data, "isosurface", params)
        if args.volumetric:
            for params in args.volumetric:
                raise ValueError("volume -volumetric not supported")

    return "cur"
parsers["volume"] = (parse_volume, parser_volume.format_usage(), parser_volume.format_help())

################################################################################

def parse_line(line, settings, state, renderer=None):
    if re.search('^\s*#', line) or len(line.strip()) == 0:
        return None

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

def parse_file(filename, settings, state=None):
    with open(filename) as fin:
        for l in fin.readlines():
            refresh = parse_line(l, settings, state)
            if state is not None:
                state.update(refresh)
