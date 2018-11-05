import argparse, numpy as np, vtk

def piecewise_linear(x, t):
    i = np.searchsorted(t[::4], x)
    if i == 0:
        return t[1:4]
    elif i == len(t)/4:
        return t[-3:]
    else:
        f = t[i*4]-x
        return f*t[(i-1)*4+1:(i-1)*4+4] + (1.0-f)*t[i*4+1:i*4+4]

config_parsers = {}

parser_colormap = argparse.ArgumentParser(prog="colormap")
parser_colormap.add_argument("-name",type=str,required=True)
parser_colormap.add_argument("colormap",nargs='+',type=float)
def parse_colormap(config, line):
    args = parser_colormap.parse_args(line.split()[1:])
    if len(args.colormap) % 4 != 0:
        raise ValueError("colormap arguments must be multiple of 4: v r g b")
    config["colormaps"][args.name] = lambda x : piecewise_linear(x, np.array(args.colormap))
config_parsers["colormap"] = parse_colormap

parser_atom_type = argparse.ArgumentParser(prog="atom_type")
parser_atom_type.add_argument("-name",type=str,required=True)
parser_atom_type.add_argument("-color",nargs=3,type=float,default=None)
parser_atom_type.add_argument("-colormap",type=str,default=None)
parser_atom_type.add_argument("-colormap_field",type=str,default=None)
parser_atom_type.add_argument("-radius",type=float,default=None)
parser_atom_type.add_argument("-radius_field",type=str,default=None)
parser_atom_type.add_argument("-opacity",type=float,default=None)
parser_atom_type.add_argument("-label_field",type=str,default=None)
def parse_atom_type(config, line):
    args = parser_atom_type.parse_args(line.split()[1:])
    if "args.name" not in config["atom_types"]:
        config["atom_types"][args.name] = {}
        config["atom_types"][args.name]["radius"] = 0.3
        config["atom_types"][args.name]["colormap_func"] = None
        config["atom_types"][args.name]["colormap_field"] = None
        config["atom_types"][args.name]["radius_field"] = None
        prop = vtk.vtkProperty()
        prop.SetOpacity(1.0)
        prop.SetSpecularColor(1.0,1.0,1.0)
        prop.SetSpecularPower(10.0)
        config["atom_types"][args.name]["prop"] = prop
    if args.color is not None:
        config["atom_types"][args.name]["prop"].SetDiffuseColor(args.color)
    if args.colormap is not None:
        config["atom_types"][args.name]["colormap_func"] = config["colormaps"][args.colormap]
    if args.colormap_field is not None:
        config["atom_types"][args.name]["colormap_field"] = args.colormap_field
    if args.opacity is not None:
        config["atom_types"][args.name]["prop"].SetOpacity(args.opacity)
    if args.radius is not None:
        config["atom_types"][args.name]["radius"] = args.radius
    if args.radius_field is not None:
        config["atom_types"][args.name]["radius_field"] = args.radius_field
    if args.label_field is not None:
        config["atom_types"][args.name]["label_field"] = args.label_field
config_parsers["atom_type"] = parse_atom_type

parser_bond_type = argparse.ArgumentParser(prog="bond_type")
parser_bond_type.add_argument("-name",type=str,required=True)
parser_bond_type.add_argument("-color",nargs=3,type=float,default=None)
parser_bond_type.add_argument("-radius",type=float,default=None)
parser_bond_type.add_argument("-opacity",type=float,default=None)
def parse_bond_type(config, line):
    args = parser_bond_type.parse_args(line.split()[1:])
    if "args.name" not in config["bond_types"]:
        config["bond_types"][args.name] = {}
        config["bond_types"][args.name]["opacity"] = 1.0
        config["bond_types"][args.name]["radius"] = 0.3
        prop = vtk.vtkProperty()
        prop.SetOpacity(1.0)
        prop.SetSpecularColor(1.0,1.0,1.0)
        prop.SetSpecularPower(10.0)
        config["bond_types"][args.name]["prop"] = prop
    if args.color is not None:
        config["bond_types"][args.name]["prop"].SetDiffuseColor(args.color)
    if args.opacity is not None:
        config["bond_types"][args.name]["prop"].SetOpacity(args.opacity)
    if args.radius is not None:
        config["bond_types"][args.name]["radius"] = args.radius
config_parsers["bond_type"] = parse_bond_type

parser_cell_box_color = argparse.ArgumentParser(prog="cell_box_color")
parser_cell_box_color.add_argument("-color",nargs=3,type=float,default=None)
def parse_bond_type(config, line):
    args = parser_cell_box_color.parse_args(line.split()[1:])
    config["cell_box_color"] = args.color
config_parsers["cell_box_color"] = parse_bond_type

parser_background_color = argparse.ArgumentParser(prog="background_color")
parser_background_color.add_argument("-color",nargs=3,type=float,default=None)
def parse_bond_type(config, line):
    args = parser_background_color.parse_args(line.split()[1:])
    config["background_color"] = args.color
config_parsers["background_color"] = parse_bond_type

def config_parse_line(config, line):
    keyword = line.split()[0]
    if keyword in config_parsers:
        config_parsers[keyword](config, line)

def config_parse_file(filename):
    config = { "atom_types" : {}, "bond_types" : {}, "colormaps" : {}, "cell_box_color" : [1.0, 1.0, 1.0], "background_color" : [0.0, 0.0, 0.0]  }
    for line in open(filename).readlines():
        config_parse_line(config, line)
    return config
