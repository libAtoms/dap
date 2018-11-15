import argparse, numpy as np, vtk
from davtk_parse_utils import ThrowingArgumentParser

def piecewise_linear(x, t):
    i = np.searchsorted(t[::4], x)
    if i == 0:
        return t[1:4]
    elif i == len(t)/4:
        return t[-3:]
    else:
        f = t[i*4]-x
        return f*t[(i-1)*4+1:(i-1)*4+4] + (1.0-f)*t[i*4+1:i*4+4]

class UnknownSettingsKeywordError(Exception):
    pass

class DavTKAtomTypes(object):
    def __init__(self):
        self.types = {}
        self.autogen_used = 0
        self.colors = []

        colorScheme = vtk.vtkColorSeries()
        for i in [1, 2, 3]:
            colorScheme.SetColorSchemeByName("Brewer Qualitative Set{}".format(i))
            for i in range(colorScheme.GetNumberOfColors()):
                c = colorScheme.GetColor(i)
                if np.any(np.array(c) != 0):
                    self.colors.append(np.array(c)/255.0)

    def __getitem__(self, key):
        if key in self.types:
            return self.types[key]
        else:
            return self.add_autogen(key)

    def add_autogen(self, name):
        print "autogenerating with color ",self.colors[self.autogen_used]
        self.set_type(name=name, color=self.colors[self.autogen_used], radius = 0.5, opacity=1.0)
        self.autogen_used += 1
        return self.types[name]

    def set_type(self, name, color=None, colormap=None, radius=None, radius_field=None, opacity=None, label=None, bonding_radius=None, colormaps=None):
        if name not in self.types:
            self.types[name] = {}
            self.types[name]["color"] = None
            self.types[name]["colormap_func"] = None
            self.types[name]["colormap_field"] = None
            self.types[name]["radius"] = None
            self.types[name]["radius_field"] = None
            self.types[name]["opacity"] = None
            self.types[name]["label"] = None
            self.types[name]["bonding_radius"] = None
            self.types[name]["prop"] = vtk.vtkProperty()
            self.types[name]["prop"].SetOpacity(1.0)
            self.types[name]["prop"].SetSpecularColor(1.0,1.0,1.0)
            self.types[name]["prop"].SetSpecularPower(10.0)
        if color is not None:
            if colormap is not None:
                raise ValueError("got color and colormap")
            self.types[name]["color"] = color
            self.types[name]["prop"].SetColor(color)
        if colormap is not None:
            if color is not None:
                raise ValueError("got color and colormap")
            self.types[name]["colormap_func"] = colormaps[colormap[0]]
            self.types[name]["colormap_field"] = colormap[1]
        if radius is not None:
            if radius_field is not None:
                raise ValueError("got radius and radius_field")
            self.types[name]["radius"] = radius
            self.types[name]["radius_field"] = None
        if radius_field is not None:
            if radius is not None:
                raise ValueError("got radius and radius_field")
            self.types[name]["radius"] = radius
            self.types[name]["radius_field"] = radius_field
        if opacity is not None:
            self.types[name]["opacity"] = opacity
            self.types[name]["prop"].SetOpacity(opacity)
        if label is not None:
            if label == "NONE" or label == "_":
                self.types[name]["label"] = None
            else:
                self.types[name]["label"] = label
        if bonding_radius is not None:
            self.types[name]["bonding_radius"] = bonding_radius

class DavTKSettings(object):
    def __init__(self):
        self.settings = { "atom_types" : DavTKAtomTypes(), "bond_types" : {}, "colormaps" : {},
            "cell_box_color" : [1.0, 1.0, 1.0], "background_color" : [0.0, 0.0, 0.0],
            "picked_color" : [1.0, 1.0, 0.0], 
            "config_n_text_color" : [1.0, 1.0, 1.0], "config_n_text_fontsize" : 36,
            "label_text_color" : [1.0, 1.0, 1.0], "label_text_fontsize" : 24 }

        self.parsers = {}

        self.parser_colormap = ThrowingArgumentParser(prog="colormap", description="args: V R G B ...")
        self.parser_colormap.add_argument("-name",type=str, required=True)
        self.parser_colormap.add_argument("colormap",nargs='*',type=float, metavar=("X"))
        self.parsers["colormap"] = self.parse_colormap

        self.parser_atom_type = ThrowingArgumentParser(prog="atom_type")
        self.parser_atom_type.add_argument("-name",type=str, required=True)
        self.parser_atom_type.add_argument("-color",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parser_atom_type.add_argument("-colormap",nargs=2,type=str,default=None, metavar=("COLORMAP","FIELD"))
        self.parser_atom_type.add_argument("-radius",type=float,default=None)
        self.parser_atom_type.add_argument("-radius_field",type=str,default=None)
        self.parser_atom_type.add_argument("-opacity",type=float,default=None)
        self.parser_atom_type.add_argument("-label",type=str,default=None)
        self.parser_atom_type.add_argument("-bonding_radius",type=float,default=None)
        self.parsers["atom_type"] = self.parse_atom_type

        self.parser_bond_type = ThrowingArgumentParser(prog="bond_type")
        self.parser_bond_type.add_argument("-name",type=str, required=True)
        self.parser_bond_type.add_argument("-color",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parser_bond_type.add_argument("-radius",type=float,default=None)
        self.parser_bond_type.add_argument("-opacity",type=float,default=None)
        self.parsers["bond_type"] = self.parse_bond_type

        self.parser_cell_box_color = ThrowingArgumentParser(prog="cell_box_color")
        self.parser_cell_box_color.add_argument("color",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parsers["cell_box_color"] = self.parse_cell_box_color

        self.parser_config_n_text = ThrowingArgumentParser(prog="config_n_text")
        self.parser_config_n_text.add_argument("-color",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parser_config_n_text.add_argument("-fontsize",type=int,default=None)
        self.parsers["config_n_text"] = self.parse_config_n_text

        self.parser_label_text = ThrowingArgumentParser(prog="label_text")
        self.parser_label_text.add_argument("-color",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parser_label_text.add_argument("-fontsize",type=int,default=None)
        self.parsers["label_text"] = self.parse_label_text

        self.parser_picked_color = ThrowingArgumentParser(prog="picked_color")
        self.parser_picked_color.add_argument("color",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parsers["picked_color"] = self.parse_picked_color

        self.parser_background_color = ThrowingArgumentParser(prog="background_color")
        self.parser_background_color.add_argument("color",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parsers["background_color"] = self.parse_background_color

        # properties
        # 3D Actor properties
        for f in ["cell_box","picked"]:
            prop = vtk.vtkProperty()
            prop.SetOpacity(1.0)
            prop.SetColor(self.settings[f+"_color"])
            self.settings[f+"_prop"] = prop

        # make picked very flat
        self.settings["picked_prop"].SetAmbient(0.6)
        self.settings["picked_prop"].SetDiffuse(0.4)

        # text properties
        for f in ["config_n", "label"]:
            prop = vtk.vtkTextProperty()
            self.settings[f+"_text_prop" ] = prop
            prop.SetOpacity(1.0)
            prop.SetColor(self.settings[f+"_text_color"])
            prop.SetFontSize(self.settings[f+"_text_fontsize"])

    def __getitem__(self,key):
        return self.settings[key]

    def parse_colormap(self, args):
        args = self.parser_colormap.parse_args(args)
        if len(args.colormap) % 4 != 0:
            raise ValueError("colormap arguments must be multiple of 4: v r g b")
        self.settings["colormaps"][args.name] = lambda x : piecewise_linear(x, np.array(args.colormap))

    def parse_atom_type(self, args):
        args = self.parser_atom_type.parse_args(args)
        self.settings["atom_types"].set_type(args.name, args.color, args.colormap, args.radius, args.radius_field, args.opacity, args.label, args.bonding_radius, self.settings["colormaps"])
        return None

    def parse_bond_type(self, args):
        refresh = None
        args = self.parser_bond_type.parse_args(args)
        if len(self.settings["bond_types"].keys()) == 0:
            self.settings["default_bond_type"] = args.name
        if args.name not in self.settings["bond_types"]:
            self.settings["bond_types"][args.name] = {}
            self.settings["bond_types"][args.name]["opacity"] = 1.0
            self.settings["bond_types"][args.name]["radius"] = 0.3
            prop = vtk.vtkProperty()
            prop.SetOpacity(1.0)
            prop.SetSpecularColor(1.0,1.0,1.0)
            prop.SetSpecularPower(10.0)
            self.settings["bond_types"][args.name]["prop"] = prop
            refresh = "all"
        if args.color is not None:
            self.settings["bond_types"][args.name]["prop"].SetColor(args.color)
        if args.opacity is not None:
            self.settings["bond_types"][args.name]["prop"].SetOpacity(args.opacity)
        if args.radius is not None:
            refresh = "all"
            self.settings["bond_types"][args.name]["radius"] = args.radius
        return refresh

    def parse_cell_box_color(self, args):
        args = self.parser_cell_box_color.parse_args(args)
        self.settings["cell_box_color"] = args.color
        self.settings["cell_box_prop"].SetColor(self.settings["cell_box_color"])
        return None

    def parse_config_n_text(self, args):
        args = self.parser_config_n_text.parse_args(args)
        if args.color is not None:
            self.settings["config_n_text_color"] = args.color
        if args.fontsize is not None:
            self.settings["config_n_text_fontsize"] = args.fontsize
        self.settings["config_n_text_prop"].SetColor(self.settings["config_n_text_color"])
        self.settings["config_n_text_prop"].SetFontSize(self.settings["config_n_text_fontsize"])
        return None

    def parse_label_text(self, args):
        args = self.parser_label_text.parse_args(args)
        if args.color is not None:
            self.settings["label_text_color"] = args.color
        if args.fontsize is not None:
            self.settings["label_text_fontsize"] = args.fontsize
        self.settings["label_text_prop"].SetColor(self.settings["label_text_color"])
        self.settings["label_text_prop"].SetFontSize(self.settings["label_text_fontsize"])
        return None

    def parse_picked_color(self, args):
        args = self.parser_picked_color.parse_args(args)
        self.settings["picked_color"] = args.color
        self.settings["picked_prop"].SetColor(self.settings["picked_color"])
        return None

    def parse_background_color(self, args):
        args = self.parser_background_color.parse_args(args)
        self.settings["background_color"] = args.color
        return "renderer"

    def parse_line(self, line):
        args = line.split()
        if args[0] in self.parsers:
            return self.parsers[args[0]](args[1:])
        else:
            raise UnknownSettingsKeywordError(args[0])
