from __future__ import print_function
from davtk.parse_utils import ThrowingArgumentParser

interactive_parsers = {}

parser_next = ThrowingArgumentParser(prog="next")
parser_next.add_argument("-n",type=int,default=1)
def parse_next(davtk_state, line):
    args = parser_next.parse_args(line.split()[1:])
    davtk_state.show_frame(dframe = args.n)
interactive_parsers["next"] = parse_next

parser_prev = ThrowingArgumentParser(prog="prev")
parser_prev.add_argument("-n",type=int,default=1)
def parse_prev(davtk_state, line):
    args = parser_prev.parse_args(line.split()[1:])
    davtk_state.show_frame(dframe = -args.n)
interactive_parsers["prev"] = parse_prev
