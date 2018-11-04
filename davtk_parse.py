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

parser_colormap = argparse.ArgumentParser(prog="colormap")
parser_colormap.add_argument("--name",type=str,required=True)
parser_colormap.add_argument("colormap",nargs='+',type=float)
def parse_colormap(defaults, line):
    args = parser_colormap.parse_args(line.split()[1:])
    if len(args.colormap) % 4 != 0:
        raise ValueError("colormap arguments must be multiple of 4: v r g b")
    defaults["colormaps"][args.name] = lambda x : piecewise_linear(x, np.array(args.colormap))

parser_atom_type = argparse.ArgumentParser(prog="atom_type")
parser_atom_type.add_argument("--name",type=str,required=True)
parser_atom_type.add_argument("--color",nargs=3,type=float,default=None)
parser_atom_type.add_argument("--colormap",type=str,default=None)
parser_atom_type.add_argument("--colormap_field",type=str,default=None)
parser_atom_type.add_argument("--radius",type=float,default=None)
parser_atom_type.add_argument("--radius_field",type=str,default=None)
parser_atom_type.add_argument("--opacity",type=float,default=None)
def parse_atom_type(defaults, line):
    args = parser_atom_type.parse_args(line.split()[1:])
    if "args.name" not in defaults["atom_types"]:
        defaults["atom_types"][args.name] = {}
        defaults["atom_types"][args.name]["radius"] = 0.3
        defaults["atom_types"][args.name]["colormap_func"] = None
        defaults["atom_types"][args.name]["colormap_field"] = None
        defaults["atom_types"][args.name]["radius_field"] = None
        prop = vtk.vtkProperty()
        prop.SetOpacity(1.0)
        prop.SetSpecularColor(1.0,1.0,1.0)
        prop.SetSpecularPower(10.0)
        defaults["atom_types"][args.name]["prop"] = prop
    if args.color is not None:
        defaults["atom_types"][args.name]["prop"].SetDiffuseColor(args.color)
    if args.colormap is not None:
        defaults["atom_types"][args.name]["colormap_func"] = defaults["colormaps"][args.colormap]
    if args.colormap_field is not None:
        defaults["atom_types"][args.name]["colormap_field"] = args.colormap_field
    if args.opacity is not None:
        defaults["atom_types"][args.name]["prop"].SetOpacity(args.opacity)
    if args.radius is not None:
        defaults["atom_types"][args.name]["radius"] = args.radius
    if args.radius_field is not None:
        defaults["atom_types"][args.name]["radius_field"] = args.radius_field

parser_bond_type = argparse.ArgumentParser(prog="bond_type")
parser_bond_type.add_argument("--name",type=str,required=True)
parser_bond_type.add_argument("--color",nargs=3,type=float,default=None)
parser_bond_type.add_argument("--radius",type=float,default=None)
parser_bond_type.add_argument("--opacity",type=float,default=None)
def parse_bond_type(defaults, line):
    args = parser_bond_type.parse_args(line.split()[1:])
    if "args.name" not in defaults["bond_types"]:
        defaults["bond_types"][args.name] = {}
        defaults["bond_types"][args.name]["opacity"] = 1.0
        defaults["bond_types"][args.name]["radius"] = 0.3
        prop = vtk.vtkProperty()
        prop.SetOpacity(1.0)
        prop.SetSpecularColor(1.0,1.0,1.0)
        prop.SetSpecularPower(10.0)
        defaults["bond_types"][args.name]["prop"] = prop
    if args.color is not None:
        defaults["bond_types"][args.name]["prop"].SetDiffuseColor(args.color)
    if args.opacity is not None:
        defaults["bond_types"][args.name]["prop"].SetOpacity(args.opacity)
    if args.radius is not None:
        defaults["bond_types"][args.name]["radius"] = args.radius

def parse_line(defaults, line):
    if line.startswith("atom_type"):
        parse_atom_type(defaults, line)
    elif line.startswith("bond_type"):
        parse_bond_type(defaults, line)
    elif line.startswith("colormap"):
        parse_colormap(defaults,line)
    elif line.startswith("cell_box_color"):
        defaults["cell_box_color"] = [float(c) for c in line.split()[1:]]
    elif line.startswith("background_color"):
        defaults["background_color"] = [float(c) for c in line.split()[1:]]

def parse_defaults(filename):
    defaults = { "atom_types" : {}, "bond_types" : {}, "colormaps" : {}, "cell_box_color" : [1.0, 1.0, 1.0], "background_color" : [0.0, 0.0, 0.0]  }
    for line in open(filename).readlines():
        parse_line(defaults, line)
    return defaults
