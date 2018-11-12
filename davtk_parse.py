from davtk_parse_utils import ThrowingArgumentParser
from davtk_settings import UnknownSettingsKeywordError

parsers = {}

parser_next = ThrowingArgumentParser(prog="next",description="go forward a number of frames")
parser_next.add_argument("-n",type=int,default=1,help="number of frames to change")
def parse_next(davtk_state, renderer, args):
    args = parser_next.parse_args(args)
    davtk_state.set_shown_frame(dframe = args.n, renderer=renderer)
    return None
parsers["next"] = (parse_next, parser_next.format_usage(), parser_next.format_help())

parser_prev = ThrowingArgumentParser(prog="prev", description="go back a number of frames")
parser_prev.add_argument("-n",type=int,default=1, help="number of frames to change")
def parse_prev(davtk_state, renderer, args):
    args = parser_prev.parse_args(args)
    davtk_state.set_shown_frame(dframe = -args.n, renderer=renderer)
    return None
parsers["prev"] = (parse_prev, parser_prev.format_usage(), parser_prev.format_help())

parser_unpick = ThrowingArgumentParser(prog="unpick")
parser_unpick.add_argument("-all",action="store_true")
def parse_unpick(davtk_state, renderer, args):
    args = parser_unpick.parse_args(args)
    if args.all:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]
    for at in ats:
        at.arrays["_vtk_picked"][:] = False
    return "all"
parsers["unpick"] = (parse_unpick, parser_unpick.format_usage(), parser_unpick.format_help())

parser_pick = ThrowingArgumentParser(prog="pick")
parser_pick.add_argument("-all",action="store_true")
parser_pick.add_argument("n",type=int,nargs='+')
def parse_pick(davtk_state, renderer, args):
    args = parser_pick.parse_args(args)
    if args.all:
        ats = davtk_state.at_list
    else:
        ats = [davtk_state.cur_at()]
    for at in ats:
        at.arrays["_vtk_picked"][args.n] = True
    return "all"
parsers["pick"] = (parse_pick, parser_pick.format_usage(), parser_pick.format_help())

parser_delete = ThrowingArgumentParser(prog="delete",description="delete atoms, picked unless indices are listed")
parser_delete.add_argument("-all",action="store_true",help="apply to all frames")
parser_delete.add_argument("n",type=int,nargs='*',help="delete by these indices ")
def parse_delete(davtk_state, renderer, args):
    args = parser_delete.parse_args(args)
    if args.all:
        frame_list=None
    else:
        frame_list="cur"
    if len(args.n) == 0:
        davtk_state.delete(atoms="picked", frames=frame_list)
    else:
        davtk_state.delete(atoms=args.n, frames=frame_list)
    return None

parsers["delete"] = (parse_delete, parser_delete.format_usage(), parser_delete.format_help())

def parse_usage(davtk_state, renderer, args):
    for keyword in sorted(parsers.keys()):
        print keyword, parsers[keyword][1],
parsers["usage"] = (parse_usage, "usage: usage\n", "usage: usage\n")

def parse_help(davtk_state, renderer, args):
    for keyword in sorted(parsers.keys()):
        print "--------------------------------------------------------------------------------"
        print keyword, parsers[keyword][2],
parsers["help"] = (parse_help, "usage: help\n", "usage: help\n")

def parse_line(line, settings, state, renderer=None):
    try:
        return settings.parse_line(line)
    except UnknownSettingsKeywordError:
        pass

    args = line.split()
    return parsers[args[0]][0](state, renderer, args[1:])

def parse_file(filename, settings, state=None):
    with open(filename) as fin:
        for l in fin.readlines():
            if l.startswith("#"):
                continue
            refresh = parse_line(l, settings, state)
            if state is not None:
                state.update(refresh)
