from __future__ import print_function
import argparse, numpy as np, vtk
import re, sys
from davtk.parse_utils import ThrowingArgumentParser

def piecewise_linear(x, t):
    i = np.searchsorted(t[::4], x)
    if i == 0:
        return t[1:4]
    elif i == len(t)/4:
        return t[-3:]
    else:
        f = t[i*4]-x
        return f*t[(i-1)*4+1:(i-1)*4+4] + (1.0-f)*t[i*4+1:i*4+4]

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
        print("autogenerating with color ",self.colors[self.autogen_used])
        self.set_type(name=name, color=self.colors[self.autogen_used], radius = 0.5, opacity=1.0)
        self.autogen_used += 1
        return self.types[name]

    def get_all(self):
        data = {}
        for name in self.types:
            t = self.types[name]
            data[name] = (t["color"],t["colormap_field"],t["radius"],t["radius_field"],
                          t["opacity"],t["label"],t["bonding_radius"])
        return data

    def set_type(self, name, color=None, colormap=None, radius=None, radius_field=None,
                 opacity=None, label=None, bonding_radius=None, colormaps=None):
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
            self.types[name]["colormap_func"] = colormaps[colormap[0]]['f']
            self.types[name]["colormap_field"] = colormap[1]
        if radius is not None:
            if radius_field is not None:
                raise ValueError("got radius and radius_field")
            self.types[name]["radius"] = radius
            self.types[name]["radius_field"] = None
        if radius_field is not None:
            if radius is not None:
                raise ValueError("got radius_field and radius")
            self.types[name]["radius"] = None
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
            "label_text_color" : [1.0, 1.0, 1.0], "label_text_fontsize" : 24,
            "frame_step" : 1, "legend" : { 'show' : False, 'position' : np.array([-100,-100]), 'spacing' : 100 }
            }

        self.parsers = {}
        
        self.parser_print_settings = ThrowingArgumentParser(prog="print_settings",description="print settings")
        self.parser_print_settings.add_argument("-keyword_regexp",type=str)
        self.parsers["print_settings"] = (self.parse_print_settings, self.parser_print_settings.format_usage(), self.parser_print_settings.format_help(), None)

        self.parser_legend = ThrowingArgumentParser(prog="legend",description="control legend, toggle by default")
        group = self.parser_legend.add_mutually_exclusive_group()
        group.add_argument("-on",action='store_true',help="enable legend")
        group.add_argument("-off",action='store_true',help="disable legend")
        group = self.parser_legend.add_mutually_exclusive_group()
        group.add_argument("-position",type=int,nargs=2,help="position relative to bottom left corner of display"+
                                                             " (negative values relative to top right)", default=None)
        group.add_argument("-offset",type=int,nargs=2,help="offset relative to current position", default=None)
        self.parser_legend.add_argument("-spacing",type=int,help="spacing between rows", default=None)
        self.parsers["legend"] = (self.parse_legend, self.parser_legend.format_usage(), self.parser_legend.format_help(), self.write_legend)

        self.parser_step = ThrowingArgumentParser(prog="step",description="number of frames to skip in +/- and prev/next")
        self.parser_step.add_argument("n",type=int,help="number of frames to step")
        self.parsers["step"] = (self.parse_step, self.parser_step.format_usage(), self.parser_step.format_help(), self.write_step)

        self.parser_colormap = ThrowingArgumentParser(prog="colormap", description="repeated sequence of groups of 4 numbers: V R G B ...")
        self.parser_colormap.add_argument("name",type=str)
        self.parser_colormap.add_argument("-P",dest="colormap", nargs=4,action='append',type=float, metavar=('V','R','G','B'))
        self.parsers["colormap"] = (self.parse_colormap, self.parser_colormap.format_usage(), self.parser_colormap.format_help(), self.write_colormap)

        self.parser_atom_type = ThrowingArgumentParser(prog="atom_type")
        self.parser_atom_type.add_argument("name",type=str)
        group = self.parser_atom_type.add_mutually_exclusive_group()
        group.add_argument("-color","-c",nargs=3,type=float,default=None, metavar=("R","G","B"))
        group.add_argument("-colormap",nargs=2,type=str,default=None, metavar=("COLORMAP","FIELD"))
        group = self.parser_atom_type.add_mutually_exclusive_group()
        group.add_argument("-radius",type=float,default=None)
        group.add_argument("-radius_field",type=str,nargs=2,metavar=("RADIUS_FIELD","FACTOR"),default=None)
        self.parser_atom_type.add_argument("-opacity",type=float,default=None)
        self.parser_atom_type.add_argument("-label_field",type=str,default=None)
        self.parser_atom_type.add_argument("-bonding_radius",type=float,default=None)
        self.parsers["atom_type"] = (self.parse_atom_type, self.parser_atom_type.format_usage(), self.parser_atom_type.format_help(), self.write_atom_type)

        self.parser_bond_type = ThrowingArgumentParser(prog="bond_type")
        self.parser_bond_type.add_argument("name",type=str)
        self.parser_bond_type.add_argument("-color","-c",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parser_bond_type.add_argument("-radius",type=float,default=None)
        self.parser_bond_type.add_argument("-opacity",type=float,default=None)
        self.parser_bond_type.add_argument("-default",action='store_true')
        self.parsers["bond_type"] = (self.parse_bond_type, self.parser_bond_type.format_usage(), self.parser_bond_type.format_help(), self.write_bond_type)

        self.parser_cell_box_color = ThrowingArgumentParser(prog="cell_box_color")
        self.parser_cell_box_color.add_argument("R",type=float)
        self.parser_cell_box_color.add_argument("G",type=float)
        self.parser_cell_box_color.add_argument("B",type=float)
        self.parsers["cell_box_color"] = (self.parse_cell_box_color, self.parser_cell_box_color.format_usage(), self.parser_cell_box_color.format_help(), self.write_cell_box_color)

        self.parser_picked_color = ThrowingArgumentParser(prog="picked_color")
        self.parser_picked_color.add_argument("R",type=float)
        self.parser_picked_color.add_argument("G",type=float)
        self.parser_picked_color.add_argument("B",type=float)
        self.parsers["picked_color"] = (self.parse_picked_color, self.parser_picked_color.format_usage(), self.parser_picked_color.format_help(), self.write_picked_color)

        self.parser_background_color = ThrowingArgumentParser(prog="background_color")
        self.parser_background_color.add_argument("R",type=float)
        self.parser_background_color.add_argument("G",type=float)
        self.parser_background_color.add_argument("B",type=float)
        self.parsers["background_color"] = (self.parse_background_color, self.parser_background_color.format_usage(), self.parser_background_color.format_help(), self.write_picked_color)

        self.parser_config_n_text = ThrowingArgumentParser(prog="config_n_text")
        self.parser_config_n_text.add_argument("-color","-c",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parser_config_n_text.add_argument("-fontsize",type=int,default=None)
        self.parsers["config_n_text"] = (self.parse_config_n_text, self.parser_config_n_text.format_usage(), self.parser_config_n_text.format_help(), self.write_config_n_text)

        self.parser_label_text = ThrowingArgumentParser(prog="label_text")
        self.parser_label_text.add_argument("-color","-c",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parser_label_text.add_argument("-fontsize",type=int,default=None)
        self.parsers["label_text"] = (self.parse_label_text, self.parser_label_text.format_usage(), self.parser_label_text.format_help(), self.write_label_text)

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

    def write(self, fout, key_re=None):
        for keyword in self.parsers:
            if (key_re is None or re.search(key_re, keyword)) and self.parsers[keyword][3] is not None:
                fout.write(self.parsers[keyword][3]())

    def parse_print_settings(self, args):
        args = self.parser_print_settings.parse_args(args)
        self.write(sys.stdout, key_re=args.keyword_regexp)
        return None

    def write_legend(self):
        args_str = 'legend '
        args_str += '-on' if self.settings["legend"]['show'] else '-off'
        args_str += ' -position {} {}'.format(self.settings["legend"]['position'][0],self.settings["legend"]['position'][1])
        args_str += ' -spacing {}'.format(self.settings["legend"]['spacing'])
        return args_str+'\n'
    def parse_legend(self, args):
        args = self.parser_legend.parse_args(args)

        if args.on:
            self.settings["legend"]['show'] = True
        elif args.off:
            self.settings["legend"]['show']= False

        if args.position is not None:
            self.settings["legend"]['position'] = np.array(args.position)
        if args.offset is not None:
            self.settings["legend"]['position'] += args.offset
        if args.spacing is not None:
            self.settings["legend"]['spacing'] = args.spacing

        if args.position is None and args.offset is None and args.spacing is None and not args.on and not args.off:
            self.settings["legend"]['show'] = not self.settings["legend"]['show']

        return "cur"

    def write_step(self):
        args_str = 'step'
        args_str += ' {}'.format(self.settings["frame_step"])
        return args_str+'\n'
    def parse_step(self, args):
        args = self.parser_step.parse_args(args)
        self.settings["frame_step"] = args.n
        return None

    def write_colormap (self):
        args_str = ''
        for colormap in self.settings["colormaps"]:
            args_str += 'colormap {}'.format(colormap)
            data = self.settings["colormaps"][colormap]["data"]
            for i in range(0,len(self.settings["colormaps"][colormap]["data"]),4):
                args_str += '   -P {} {} {} {}'.format(data[i],data[i+1],data[i+2],data[i+3])
            args_str += '\n'
        return args_str
    def parse_colormap(self, args):
        args = self.parser_colormap.parse_args(args)
        args.colormap = [item for sublist in args.colormap for item in sublist]
        if len(args.colormap) % 4 != 0:
            raise ValueError("colormap arguments must be multiple of 4: v r g b")
        self.settings["colormaps"][args.name] = { 'f' : lambda x : piecewise_linear(x, np.array(args.colormap)),
                                                  'data' : args.colormap }
        return "cur"

    def write_atom_type (self):
        args_str = ''
        all_data = self.settings["atom_types"].get_all()
        for atom_type in all_data:
            args_str += 'atom_type {}'.format(atom_type)
            (color, colormap, radius, radius_field, opacity, label_field, bonding_radius) = all_data[atom_type]
            if color is not None:
                args_str += ' -color {} {} {}'.format(color[0], color[1], color[2])
            if colormap is not None:
                args_str += ' -colormap {}'.format(colormap)
            if radius is not None:
                args_str += ' -radius {}'.format(radius)
            if radius_field is not None:
                args_str += ' -radius_field {} {}'.format(radius_field[0], radius_field[1])
            if opacity is not None:
                args_str += ' -opacity {}'.format(opacity)
            if label_field is not None:
                args_str += ' -label_field {}'.format(label_field)
            if bonding_radius is not None:
                args_str += ' -bonding_radius {}'.format(bonding_radius)
            args_str += '\n'
        return args_str
    def parse_atom_type(self, args):
        args = self.parser_atom_type.parse_args(args)
        if args.radius_field is not None:
            args.radius_field[1] = float(args.radius_field[1])
        self.settings["atom_types"].set_type(args.name, args.color, args.colormap, args.radius, args.radius_field, args.opacity, args.label_field, args.bonding_radius, self.settings["colormaps"])
        return None

    def write_bond_type (self):
        args_str = ''
        for bond_type, data in self.settings["bond_types"].items():
            args_str += 'bond_type {}'.format(bond_type)
            c = data["prop"].GetColor()
            args_str += ' -color {} {} {}'.format(c[0],c[1],c[2])
            if data["opacity"] is not None:
                args_str += ' -opacity {}'.format(data["opacity"])
            if data["radius"] is not None:
                args_str += ' -radius {}'.format(data["radius"])
            args_str += '\n'
        return args_str
    def parse_bond_type(self, args):
        refresh = None
        args = self.parser_bond_type.parse_args(args)

        if len(self.settings["bond_types"].keys()) == 0: # first bond type is initial default
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

        if args.default:
            self.settings["default_bond_type"] = args.name

        return refresh

    def write_cell_box_color(self):
        args = 'cell_box_color {} {} {}'.format(self.settings["cell_box_color"][0],
                                                self.settings["cell_box_color"][1],
                                                self.settings["cell_box_color"][2])
        return args+'\n'
    def parse_cell_box_color(self, args):
        args = self.parser_cell_box_color.parse_args(args)
        self.settings["cell_box_color"] = (args.R, args.G, args.B)
        self.settings["cell_box_prop"].SetColor(self.settings["cell_box_color"])
        return None

    def write_config_n_text(self):
        args = 'config_n_text'
        if self.settings["config_n_text_color"] is not None:
            args += ' -color {} {} {}'.format(self.settings["config_n_text_color"][0],
                                              self.settings["config_n_text_color"][1],
                                              self.settings["config_n_text_color"][2])
        if self.settings["config_n_text_fontsize"] is not None:
            args += ' -fontsize {}'.format(self.settings["config_n_text_fontsize"])
        return args+'\n'
    def parse_config_n_text(self, args):
        args = self.parser_config_n_text.parse_args(args)
        if args.color is not None:
            self.settings["config_n_text_color"] = args.color
        if args.fontsize is not None:
            self.settings["config_n_text_fontsize"] = args.fontsize
        self.settings["config_n_text_prop"].SetColor(self.settings["config_n_text_color"])
        self.settings["config_n_text_prop"].SetFontSize(self.settings["config_n_text_fontsize"])
        return None

    def write_label_text(self):
        args = 'label_text'
        if self.settings["label_text_color"] is not None:
            args += ' -color {} {} {}'.format(self.settings["label_text_color"][0],
                                              self.settings["label_text_color"][1],
                                              self.settings["label_text_color"][2])
        if self.settings["label_text_fontsize"] is not None:
            args += ' -fontsize {}'.format(self.settings["label_text_fontsize"])
        return args+'\n'
    def parse_label_text(self, args):
        args = self.parser_label_text.parse_args(args)
        if args.color is not None:
            self.settings["label_text_color"] = args.color
        if args.fontsize is not None:
            self.settings["label_text_fontsize"] = args.fontsize
        self.settings["label_text_prop"].SetColor(self.settings["label_text_color"])
        self.settings["label_text_prop"].SetFontSize(self.settings["label_text_fontsize"])
        return None

    def write_picked_color(self):
        args = 'picked_color {} {} {}'.format(self.settings["picked_color"][0],
                                              self.settings["picked_color"][1],
                                              self.settings["picked_color"][2])
        return args+'\n'
    def parse_picked_color(self, args):
        args = self.parser_picked_color.parse_args(args)
        self.settings["picked_color"] = (args.R, args.G, args.B)
        self.settings["picked_prop"].SetColor(self.settings["picked_color"])
        return None

    def write_background_color(self):
        args = 'background_color {} {} {}'.format(self.settings["background_color"][0],
                                                  self.settings["background_color"][1],
                                                  self.settings["background_color"][2])
        return args+'\n'
    def parse_background_color(self, args):
        args = self.parser_background_color.parse_args(args)
        self.settings["background_color"] = (args.R, args.G, args.B)
        return "cur"
