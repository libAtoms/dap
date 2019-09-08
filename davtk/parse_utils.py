from __future__ import print_function
import argparse

# subclass ArgumentParser to throw errors instead of exiting
class ArgumentParserError(Exception):
    pass
class ArgumentParserHelp(Exception):
    pass

class ThrowingArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ArgumentParserError(message)
    def exit(self):
        raise ArgumentParserHelp("help")

def material_dict_from_args(args):
    d = {}
    for f in ["opacity", "specular", "specular_radius", "ambient"]:
        if hasattr(args,f) and getattr(args,f) is not None:
            d[f] = getattr(args,f)
    return d

def write_material_args(material_args):
    string = ''
    string += ' -opacity {}'.format(material_args["opacity"])
    string += ' -specular {}'.format(material_args["specular"])
    string += ' -specular_radius {}'.format(material_args["specular_radius"])
    string += ' -ambient {}'.format(material_args["ambient"])
    return string

def add_material_args_to_parser(arg_parser):
    arg_parser.add_argument("-opacity",type=float,default=None)
    arg_parser.add_argument("-specular",type=float,default=None)
    arg_parser.add_argument("-specular_radius",type=float,default=None)
    arg_parser.add_argument("-ambient",type=float,default=None)

