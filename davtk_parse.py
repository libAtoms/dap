from davtk_parse_utils import ThrowingArgumentParser
from davtk_config import UnknownConfigKeywordError

parsers = {}

parser_next = ThrowingArgumentParser(prog="next")
parser_next.add_argument("-n",type=int,default=1)
def parse_next(davtk_state, renderer, args):
    args = parser_next.parse_args(args)
    davtk_state.set_shown_frame(dframe = args.n, renderer=renderer)
    return None
parsers["next"] = parse_next

parser_prev = ThrowingArgumentParser(prog="prev")
parser_prev.add_argument("-n",type=int,default=1)
def parse_prev(davtk_state, renderer, args):
    args = parser_prev.parse_args(args)
    davtk_state.set_shown_frame(dframe = -args.n, renderer=renderer)
    return None
parsers["prev"] = parse_prev

def parse_line(line, config, state, renderer=None):
    try:
        return config.parse_line(line)
    except UnknownConfigKeywordError:
        pass

    args = line.split()
    return parsers[args[0]](state, renderer, args[1:])

def parse_file(filename, config, state=None):
    with open(filename) as fin:
        for l in fin.readlines():
            if l.startswith("#"):
                continue
            refresh = parse_line(l, config, state)
            if state is not None:
                state.update(refresh)
