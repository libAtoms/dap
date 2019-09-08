from __future__ import print_function

import sys, re, os, types
import numpy as np
import ase.io
from ase.calculators.vasp import VaspChargeDensity
import ffmpeg

from davtk.parse_utils import ThrowingArgumentParser, add_material_args_to_parser
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
from davtk.vtk_utils import new_prop, update_prop

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
parser_exit = ThrowingArgumentParser(prog="exit",description="end program")
def parse_exit(davtk_state, renderer, args):
    return "exit"
parsers["exit"] = (parse_exit, parser_exit.format_usage(), parser_exit.format_help())

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
        davtk_state.update(str(frame_i))
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
            # davtk_state.update(str(frame_i))
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
    davtk_state.update(str(args.n))
    return None
parsers["go"] = (parse_go, parser_go.format_usage(), parser_go.format_help())

parser_next = ThrowingArgumentParser(prog="next",description="go forward a number of frames")
parser_next.add_argument("n",type=int,nargs='?',default=0,help="number of frames to change (default set by 'step' command)")
def parse_next(davtk_state, renderer, args):
    args = parser_next.parse_args(args)
    davtk_state.update("+{}".format(args.n if args.n > 0 else davtk_state.settings["frame_step"]))
    return None
parsers["next"] = (parse_next, parser_next.format_usage(), parser_next.format_help())

parser_prev = ThrowingArgumentParser(prog="prev", description="go back a number of frames")
parser_prev.add_argument("n",type=int,nargs='?',default=0, help="number of frames to change (default set by 'step' command)")
def parse_prev(davtk_state, renderer, args):
    args = parser_prev.parse_args(args)
    davtk_state.update("-{}".format(args.n if args.n > 0 else davtk_state.settings["frame_step"]))
    return None
parsers["prev"] = (parse_prev, parser_prev.format_usage(), parser_prev.format_help())

parser_unpick = ThrowingArgumentParser(prog="unpick", description="unpick picked atoms and/or bonds (default both)")
parser_unpick.add_argument("-all_frames", action="store_true")
parser_unpick.add_argument("-atoms", action="store_true")
parser_unpick.add_argument("-bonds", action="store_true")
def parse_unpick(davtk_state, renderer, args):
    args = parser_unpick.parse_args(args)
    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]
    none_selected = not any([args.atoms, args.bonds])
    for at in ats:
        if none_selected or args.atoms:
            if "_vtk_picked" not in at.arrays:
                continue
            at.arrays["_vtk_picked"][:] = False
        if none_selected or args.bonds:
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
        if "_vtk_picked" not in at.arrays:
            at.new_arrays("_vtk_picked",np.array([False]*len(at)))
        at.arrays["_vtk_picked"][args.n] = True
    return "cur"
parsers["pick"] = (parse_pick, parser_pick.format_usage(), parser_pick.format_help())

parser_delete = ThrowingArgumentParser(prog="delete",description="delete objects (picked by default)")
parser_delete.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_delete.add_argument("-atoms",type=int,nargs='+',help="delete by these indices ")
def parse_delete(davtk_state, renderer, args):
    args = parser_delete.parse_args(args)
    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]

    for at in ats:
        if args.atoms is not None:
            davtk_state.delete(at, atoms=args.atoms, frames=frame_list)
        else:
            davtk_state.delete(at, atoms="picked", bonds="picked", frames=frame_list)

    return "cur"
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

    davtk_state.update()

    return None
parsers["images"] = (parse_images, parser_images.format_usage(), parser_images.format_help())

parser_vectors = ThrowingArgumentParser(prog="vectors",description="Draw vectors")
parser_vectors.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_vectors.add_argument("-field",type=str,help="atom field to use for vectors (scalar or 3-vector)", default=None)
parser_vectors.add_argument("-color",type=str,nargs='+',help="R G B, or string atom (by atom) or string sign (by sign of scalars only)", default=None)
parser_vectors.add_argument("-sign_colors",type=float,nargs=6,metavar=('RU','GU','BU','RD','GD','BD'), help="colors for -name sign", default=None)
parser_vectors.add_argument("-radius",type=float,help="cylinder radius", default=None)
parser_vectors.add_argument("-scale",type=float,help="scaling factor from field value to cylinder length", default=None)
parser_vectors.add_argument("-delete",action='store_true',help="disable")
def parse_vectors(davtk_state, renderer, args):
    args = parser_vectors.parse_args(args)

    if args.color is not None:
        if len(args.color) == 1:
            if args.color[0] != "atom" and args.color[0] != "sign":
                raise ValueError("Got single word -color argument '{}', not 'atom' or 'sign'".format(args.color))
            args.color = args.color[0]
        elif len(args.color) == 3:
            try:
                args.color = [float(v) for v in args.color]
            except:
                raise ValueError("Got 3-word args.color '{}' not convertible to float".format(args.color))
        else:
            raise ValueError("Got args.color '{}' length not 1 or 3".format(args.color))

    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]

    for at in ats:
        if args.delete:
            del at.info["_vtk_vectors"]
        else:
            if "_vtk_vectors" not in at.info:
                at.info["_vtk_vectors"] = { "field" : "magmoms", "color" : "atom", "radius" : 0.1, "scale" : 1.0, "sign_colors" : [1.0, 0.0, 0.0,    0.0, 0.0, 1.0] }
            for p in ["field", "color", "sign_colors", "radius", "scale" ]:
                if getattr(args,p) is not None:
                    at.info["_vtk_vectors"][p] = getattr(args, p)
    return None
parsers["vectors"] = (parse_vectors, parser_vectors.format_usage(), parser_vectors.format_help())

parser_bond = ThrowingArgumentParser(prog="bond",description="Create bonds")
parser_bond.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_bond.add_argument("-name",type=str,help="name of bond set", default="default")
parser_bond.add_argument("-T",type=str,help="Restrict one member to given atom_type", default="*")
parser_bond.add_argument("-T2",type=str,help="Restrict other member to given atom_type", default="*")
grp = parser_bond.add_mutually_exclusive_group()
grp.add_argument("-picked", action='store_true', help="bond a pair of picked atoms with MIC")
grp.add_argument("-n", action='store', type=int, nargs=2, help="bond a pair of atoms with given IDs", default=None)
grp.add_argument("-cutoff", "-c", action='store', type=float, nargs='+', help="bond atoms with max (one argument) or min-max (two arguments) cutoffs", default=None)
grp.add_argument("-auto", action='store_true', help="bond by previously set per-atom bonding radius")
parser_bond.add_argument("-radius", type=float, help="radius of bond cylinders", default=None)
parser_bond.add_argument("-color", nargs=3, type=float, metavar=("R","G","B"), help="color for bonds", default=None)
add_material_args_to_parser(parser_bond)
group = parser_bond.add_mutually_exclusive_group()
group.add_argument("-delete", action='store_true', help="delete existing bond")
group.add_argument("-list", action='store_true', help="list existing bond")
def parse_bond(davtk_state, renderer, args):
    args = parser_bond.parse_args(args)

    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]

    creating_bonds = args.picked or args.n is not None or args.cutoff is not None or args.auto

    if creating_bonds:
        # check for conflicting args
        if args.delete or args.list:
            raise RuntimeError("can't create bonds and also -delete or -list")
        if args.cutoff is not None and len(args.cutoff) > 2:
            raise ValueError("bond got -cutoff with {} > 2 values".format(len(args.cutoff)))

    # create new property if needed (even if not creating bonds)
    if args.name not in davtk_state.bond_prop:
        prop = new_prop( types.SimpleNamespace( color = (1.0, 1.0, 1.0), opacity = 1.0,
            specular = 0.0, specular_radius = 0.1, ambient = 0.2 ) )
        prop.radius = 0.2
        davtk_state.bond_prop[args.name] = prop

    if not args.delete and not args.list:
        # update property
        update_prop(davtk_state.bond_prop[args.name], args)
        if args.radius is not None:
            davtk_state.bond_prop[args.name].radius = args.radius

    if creating_bonds:
        # actually create bonds
        for at in ats:
            if args.picked:
                if args.T != '*' or args.T2 != '*':
                    raise RuntimeError("bond by picked can't specify -T or -T2")
                davtk_state.bond(at, args.name, args.T, args.T2, "picked")
            elif args.n:
                if args.T != '*' or args.T2 != '*':
                    raise RuntimeError("bond by n (indices) can't specify -T or -T2")
                davtk_state.bond(at, args.name, args.T, args.T2, ("n", args.n))
            elif args.cutoff:
                davtk_state.bond(at, args.name, args.T, args.T2, ("cutoff", args.cutoff))
            elif args.auto:
                davtk_state.bond(at, args.name, args.T, args.T2, ("cutoff", None))
    else: # not creating, either delete or list or just modifying an existing name
        if args.delete:
            davtk_state.delete(atoms=None, bonds=args.name)
        elif args.list:
            if len(davtk_state.bond_prop) > 0:
                print("bond names: ",list(davtk_state.bond_prop.keys()))
            else:
                print("no bond names defined")
        # else: just modifying existing proprerty

    # TODO
    # look for unused bond names and remove their props

    if creating_bonds or args.delete:
        return "cur"
    elif args.radius is None:
        return "color_only"
    else:
        return "settings"
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
    return "cur"
parsers["X"] = (parse_X, parser_X.format_usage(), parser_X.format_help())

parser_read = ThrowingArgumentParser(prog="read",description="read commands from file(s)")
parser_read.add_argument("filename",type=str,nargs='+',help="filenames")
def parse_read(davtk_state, renderer, args):
    args = parser_read.parse_args(args)
    for f in args.filename:
        parse_file(f, davtk_state.settings, davtk_state)
parsers["read"] = (parse_read, parser_read.format_usage(), parser_read.format_help())

parser_override_frame_label = ThrowingArgumentParser(prog="override_frame_label",description="set per-frame label")
parser_override_frame_label.add_argument("-all_frames",action="store_true",help="apply to all frames")
parser_override_frame_label.add_argument("string",action="store",help="per-configuration string, substituting ${FIELD} with fields "+
                                                                      "in atoms.info (or 'config_n'), or _NONE_")
def parse_override_frame_label(davtk_state, renderer, args):
    args = parser_override_frame_label.parse_args(args)
    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]

    for at in ats:
        at.info["_vtk_frame_label_string"] = args.string

    return "cur"
parsers["override_frame_label"] = (parse_override_frame_label, parser_override_frame_label.format_usage(), parser_override_frame_label.format_help())

parser_override_atom_label = ThrowingArgumentParser(prog="override_atom_label",description="override atom label on a per-atom basis")
parser_override_atom_label.add_argument("-all_frames",action="store_true",help="apply to all frames")
group = parser_override_atom_label.add_mutually_exclusive_group()
group.add_argument("-picked",action='store_true',help="label picked atoms with STRING, default", default=None)
group.add_argument("-i",type=str,action='store',help="label specified atoms with STRING, comma separated set of ranges", default=None)
group.add_argument("-eval",type=str,nargs='+',action='store',help="label specified atoms with STRING, expression that evaluates "+
                                                                  "to list of indices or logical array with dim len(atoms)", default=None)
parser_override_atom_label.add_argument("string",type=str,nargs='+',action='store',help="string to set as per-atom label", default=[])
def parse_override_atom_label(davtk_state, renderer, args):
    args = parser_override_atom_label.parse_args(args)
    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]

    for at in ats:
        new_label = " ".join(args.string)
        if "_vtk_label" not in at.arrays or len(new_label) > len(at.arrays["_vtk_label"][0]):
            a = np.array([" "*len(new_label)]*len(at))
            if "_vtk_label" in at.arrays:
                a[:] = at.arrays["_vtk_label"]
                del at.arrays["_vtk_label"]
            at.new_array("_vtk_label", a)
        if args.i:
            for range_item in args.i.split(","):
                exec('at.arrays["_vtk_label"]['+range_item+'] = new_label')
        elif args.eval:
            atoms = at
            expr = " ".join(args.eval)
            exec('atoms.arrays["_vtk_label"]['+expr+'] = new_label')
        else: # -picked, or nothing selected (default)
            if "_vtk_picked" not in at.arrays:
                return None
            for i_at in range(len(at)):
                if at.arrays["_vtk_picked"][i_at]:
                    at.arrays["_vtk_label"][i_at] = new_label
    return "cur"
parsers["override_atom_label"] = (parse_override_atom_label, parser_override_atom_label.format_usage(), parser_override_atom_label.format_help())

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

# Logic of which arguments can coexist is way too messy here. Not sure how to fix.
parser_polyhedra = ThrowingArgumentParser(prog="polyhedra",description="draw coordination polyhedra")
parser_polyhedra.add_argument("-all_frames", action="store_true", help="apply to all frames")
parser_polyhedra.add_argument("-name", help="name of polyhedron set, for later reference with -delete", default=None)
parser_polyhedra.add_argument("-T", type=str, help="atom_type for polyhedron center, required",  default=None)
parser_polyhedra.add_argument("-Tn", type=str, help="atom_type for polyhedron neighbors",  default=None)
parser_polyhedra.add_argument("-color", nargs=3, type=float, metavar=("R","G","B"), help="color for polyhedra", default=None)
add_material_args_to_parser(parser_polyhedra)
group = parser_polyhedra.add_mutually_exclusive_group()
group.add_argument("-cutoff", type=float, help="center-neighbor cutoff distance",  default=None)
group.add_argument("-bond_name", help="name of bond to use for neighbors",  default=None)
group.add_argument("-delete", action='store', metavar="NAME", help="name of polyhedra to delete")
group.add_argument("-list", action='store_true', help="list existing polyhedra")
def parse_polyhedra(davtk_state, renderer, args):
    args = parser_polyhedra.parse_args(args)

    if args.all_frames:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]

    creating_polyhedra = args.cutoff is not None or args.bond_name is not None

    if args.name is None and not args.delete and not args.list:
        raise RuntimeError("-name required when creating or modifying polyhedra")

    if creating_polyhedra:
        # check for conflicting args
        if args.delete is not None or args.list:
            raise RuntimeError("can't create polyhedra and also -delete or -list")
        if args.T is None:
            raise RuntimeError("polyhedra requires -T to create new polyhedra")

    # create new property if needed
    if args.name not in davtk_state.polyhedra_prop:
        prop = new_prop( types.SimpleNamespace( color = (0.5, 0.5, 1.0), opacity = 0.5,
            specular = 0.7, specular_radius = 0.1, ambient = 0.1 ) )
        davtk_state.polyhedra_prop[args.name] = prop

    if not args.delete is not None and not args.list:
        # update property
        update_prop(davtk_state.polyhedra_prop[args.name], args)

    if creating_polyhedra:
        # actually create
        for at in ats:
            davtk_state.polyhedra(at, args.name, args.T, args.Tn, args.cutoff, args.bond_name)
    else: # not creating, either delete or list or just modifying an existing name
        if args.T is not None or args.Tn is not None:
            raise RuntimeError("polyhedra got -T but no other polyhedra creation args")
        if args.delete is not None:
            if args.name is not None:
                raise RuntimeError("polyhedra got -delete and -name")
            for at in ats:
                if "_vtk_polyhedra_"+args.delete in at.arrays:
                    del at.arrays["_vtk_polyhedra_"+args.delete]
            # TODO: look for unused polyhedra names and remove their props ?
        elif args.list:
            if args.name is not None:
                raise RuntimeError("polyhedra got -list and -name")
            if len(davtk_state.polyhedra_prop) > 0:
                names = list(set([n.replace("_vtk_polyhedra_","") for at in ats for n in at.arrays if n.startswith("_vtk_polyhedra")]))
                print("polyhedra sets:",names)
            else:
                print("no polyhedra set names defined")
        # else: just modifying existing property


    if creating_polyhedra or args.delete:
        return "cur"
    else:
        return "color_only"
parsers["polyhedra"] = (parse_polyhedra, parser_polyhedra.format_usage(), parser_polyhedra.format_help())

parser_volume = ThrowingArgumentParser(prog="volume",description="read volumetric data from file")
parser_volume.add_argument("filename",nargs='?',help="file to read from")
parser_volume.add_argument("-name",type=str, help="name to assign", default=None)
group = parser_volume.add_mutually_exclusive_group()
group.add_argument("-delete", action='store', metavar="NAME", help="name of volume to delete", default=None)
group.add_argument("-list", action='store_true', help="list existing volume representations")
group.add_argument("-isosurface",type=float,action='store',metavar="THRESHOLD", help="isosurface threshold")
group.add_argument("-volumetric",type=float,action='store',metavar="SCALE", help="volumetric value_to_opacity_factor and color")
parser_volume.add_argument("-color", nargs=3, type=float, metavar=("R","G","B"), help="color for volumetric representation", default=None)
add_material_args_to_parser(parser_volume)
def parse_volume(davtk_state, renderer, args):
    args_list = args
    args = parser_volume.parse_args(args)

    if args.delete is not None:
        davtk_state.delete_volume_rep(davtk_state.cur_at(), args.delete)
    elif args.list:
        print("Defined volume names: "+" ".join(davtk_state.cur_at().volume_reps.keys()))
    else:
        creating_rep = args.isosurface is not None or args.volumetric is not None
        if args.filename is None and creating_rep:
            raise RuntimeError("volume got no filename")
        if args.name is None:
            if creating_rep:
                args.name = args.filename
            else:
                raise RuntimeError("not creating, deleting, or listing, need -name to modify")

        # create new property if needed
        if args.name not in davtk_state.volume_rep_prop:
            prop = new_prop( types.SimpleNamespace( color = (0.5, 0.5, 1.0), opacity = 0.5,
                specular = 0.7, specular_radius = 0.1, ambient = 0.1 ) )
            davtk_state.volume_rep_prop[args.name] = prop

        update_prop(davtk_state.volume_rep_prop[args.name], args)

        if creating_rep:
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

            if args.isosurface is not None:
                davtk_state.add_volume_rep(args.name, data, "isosurface", (args.isosurface,), ["volume"] + args_list)
            if args.volumetric is not None:
                raise RuntimeError("volume -volumetric not supported")
        # else: just modifying property

    return "cur"
parsers["volume"] = (parse_volume, parser_volume.format_usage(), parser_volume.format_help())

parser_view = ThrowingArgumentParser(prog="view",description="set view position and orientation")
parser_view.add_argument("-lattice",action='store_true',help="use lattice A1 A2 A3 instead of cartesian X Y Z directions")
parser_view.add_argument("-dir",nargs=6,type=float,metavar=("ALONG_X","ALONG_Y","ALONG_Z","UP_X","UP_Y","UP_Z"), help="view and up directions", default=None)
parser_view.add_argument("-mag",type=float,help="view magnification (relative to current)")
def parse_view(davtk_state, renderer, args):
    args = parser_view.parse_args(args)

    if args.mag <= 0:
        raise ValueError("Can't set magnification <= 0")

    davtk_state.set_view(args.dir, args.lattice, args.mag)

    return "settings"
parsers["view"] = (parse_view, parser_view.format_usage(), parser_view.format_help())
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
