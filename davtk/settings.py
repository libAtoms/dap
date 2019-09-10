from __future__ import print_function
import argparse, numpy as np, vtk
import re, sys
from davtk.parse_utils import ThrowingArgumentParser, write_material_args, add_material_args_to_parser
from davtk.vtk_utils import update_prop
import types

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
        print("autogenerating {} with color {}".format(name, self.colors[self.autogen_used]))
        self.set_type(name=name, color=self.colors[self.autogen_used], radius = 0.5)
        self.autogen_used += 1
        return self.types[name]

    def get_all(self):
        data = {}
        for name in self.types:
            t = self.types[name]
            data[name] = (t["color"],t["colormap"],t["radius"],t["radius_field"],
                          {k : t[k] for k in ["opacity","specular","specular_radius","ambient"]},
                          t["bonding_radius"])
        return data

    def set_type(self, name, color=None, colormap=None, radius=None, radius_field=None,
                 specular=None, specular_radius=None, ambient=None, opacity=None,
                 bonding_radius=None):
        if name not in self.types:
            self.types[name] = {}
            self.types[name]["color"] = None
            self.types[name]["colormap"] = None
            self.types[name]["radius"] = None
            self.types[name]["radius_field"] = None
            self.types[name]["opacity"] = 1.0
            self.types[name]["specular"] = 0.8
            self.types[name]["specular_radius"] = 0.1
            self.types[name]["ambient"] = 0.3
            self.types[name]["bonding_radius"] = None
            self.types[name]["prop"] = vtk.vtkProperty()
        if color is not None:
            if colormap is not None:
                raise ValueError("got color and colormap")
            self.types[name]["color"] = color
            self.types[name]["colormap"] = None
        if specular is not None:
            self.types[name]["specular"] = specular
        if specular_radius is not None:
            self.types[name]["specular_radius"] = specular_radius
        if ambient is not None:
            self.types[name]["ambient"] = ambient
        if colormap is not None:
            if color is not None:
                raise ValueError("got color and colormap")
            self.types[name]["colormap"] = (colormap[0], colormap[1])
            self.types[name]["color"] = None
        if radius is not None:
            if radius_field is not None:
                raise ValueError("got radius and radius_field")
            self.types[name]["radius"] = radius
            self.types[name]["radius_field"] = None
        if radius_field is not None:
            if radius is not None:
                raise ValueError("got radius_field and radius")
            self.types[name]["radius_field"] = radius_field
            self.types[name]["radius"] = None
        if opacity is not None:
            self.types[name]["opacity"] = opacity
        if bonding_radius is not None:
            self.types[name]["bonding_radius"] = bonding_radius

        update_prop(self.types[name]["prop"], types.SimpleNamespace( color=color, specular=specular, 
                    specular_radius = specular_radius, ambient=ambient, opacity=opacity ) )

class DavTKSettings(object):
    def __init__(self):
        self.data = { "atom_types" : DavTKAtomTypes(), "colormaps" : {},
            "cell_box" : { "color" : [1.0, 1.0, 1.0], "opacity" : 1.0, "line_width" : 2.0, "prop" : None }, "background_color" : [0.0, 0.0, 0.0],
            "picked" : { "color" : [1.0, 1.0, 0.0], "opacity" : 1.0,  "prop" : None },
            "frame_label" : { "color" : [1.0, 1.0, 1.0], "fontsize" : 36, "prop" : None, "string" : "${config_n}", "show" : True },
            "atom_label" : { "color" : [1.0, 1.0, 1.0], "fontsize" : 24 , "prop" : None, "string" : "$${ID}", "show" : False},
            "frame_step" : 1, "legend" : { 'show' : False, 'position' : np.array([-10,-10]), 'spacing' : 1.0, 'sphere_scale' : 1.0 },
            'atom_type_field' : 'Z'
            }

        self.parsers = {}

        self.parser_atom_type_field = ThrowingArgumentParser(prog="atom_type_field",description="ASE at.arrays field to use for atom type")
        self.parser_atom_type_field.add_argument("field",type=str,help="name of field ('Z' for atomic numbers, 'species' for chemical symbols", default='Z')
        self.parsers["atom_type_field"] = (self.parse_atom_type_field, self.parser_atom_type_field.format_usage(),
                                           self.parser_atom_type_field.format_help(), self.write_atom_type_field)

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
        self.parser_legend.add_argument("-spacing",type=float,help="multiplier for spacing between rows", default=None)
        self.parser_legend.add_argument("-sphere_scale",action='store',type=float,help="scaling factor for sphere radius", default=None)
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
        group.add_argument("-radius","-rad","-r",type=float,default=None)
        group.add_argument("-radius_field",type=str,nargs=2,metavar=("RADIUS_FIELD","FACTOR"),default=None)
        self.parser_atom_type.add_argument("-bonding_radius",type=float,default=None)
        add_material_args_to_parser(self.parser_atom_type)
        self.parsers["atom_type"] = (self.parse_atom_type, self.parser_atom_type.format_usage(), self.parser_atom_type.format_help(), self.write_atom_type)

        self.parser_cell_box = ThrowingArgumentParser(prog="cell_box")
        self.parser_cell_box.add_argument("-color",nargs=3,type=float,metavar=['R','G','B'], default=None)
        self.parser_cell_box.add_argument("-opacity",type=float,default=None)
        self.parser_cell_box.add_argument("-width",type=float,default=None)
        self.parsers["cell_box"] = (self.parse_cell_box, self.parser_cell_box.format_usage(), self.parser_cell_box.format_help(),
                                          self.write_cell_box)

        self.parser_picked = ThrowingArgumentParser(prog="picked")
        self.parser_picked.add_argument("-color",nargs=3,type=float,metavar=['R','G','B'])
        self.parsers["picked"] = (self.parse_picked, self.parser_picked.format_usage(), self.parser_picked.format_help(),
                                        self.write_picked)

        self.parser_background_color = ThrowingArgumentParser(prog="background_color")
        self.parser_background_color.add_argument("color",nargs=3,type=float,metavar=['R','G','B'])
        self.parsers["background_color"] = (self.parse_background_color, self.parser_background_color.format_usage(), self.parser_background_color.format_help(),
                                            self.write_background_color)

        self.parser_frame_label = ThrowingArgumentParser(prog="frame_label")
        self.parser_frame_label.add_argument("-string","-s", type=str,nargs='+',
            help="string, evaluating $( EXPRESSION ) and substituting ${STRING} with fields in atoms.info (or 'config_n'), or _NONE_", default=None)
        self.parser_frame_label.add_argument("-color","-c", nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parser_frame_label.add_argument("-fontsize",type=int,default=None)
        group = self.parser_frame_label.add_mutually_exclusive_group()
        group.add_argument("-on", action='store_true')
        group.add_argument("-off", action='store_true')
        self.parsers["frame_label"] = (self.parse_frame_label, self.parser_frame_label.format_usage(), self.parser_frame_label.format_help(),
                                       self.write_frame_label)

        self.parser_atom_label = ThrowingArgumentParser(prog="atom_label")
        self.parser_atom_label.add_argument("-string",type=str,help="string to use for label, evaluating $( EXPRESSION ) and substituting $${STRING} with fields in atoms.arrays "+
                                                                   "(or 'ID' for number of atom, 'Z' for atomic number, 'species' for chemical symbol), "+
                                                                   "or '_NONE_'", default=None)
        self.parser_atom_label.add_argument("-color","-c",nargs=3,type=float,default=None, metavar=("R","G","B"))
        self.parser_atom_label.add_argument("-fontsize",type=int,default=None)
        group = self.parser_atom_label.add_mutually_exclusive_group()
        group.add_argument("-on", action='store_true')
        group.add_argument("-off", action='store_true')
        self.parsers["atom_label"] = (self.parse_atom_label, self.parser_atom_label.format_usage(), self.parser_atom_label.format_help(),
                                      self.write_atom_label)

        # properties
        # 3D Actor properties
        for f in ["cell_box","picked"]:
            prop = vtk.vtkProperty()
            prop.SetColor(self.data[f]["color"])
            prop.SetOpacity(self.data[f]["opacity"])
            if "line_width" in self.data[f]:
                prop.SetLineWidth(self.data[f]["line_width"])
            self.data[f]["prop"] = prop

        # make picked very flat
        self.data["picked"]["prop"].SetAmbient(0.6)
        self.data["picked"]["prop"].SetDiffuse(0.4)

        # text properties
        for f in ["frame_label", "atom_label"]:
            prop = vtk.vtkTextProperty()
            self.data[f]["prop"] = prop
            prop.SetOpacity(1.0)
            prop.SetColor(self.data[f]["color"])
            prop.SetFontSize(self.data[f]["fontsize"])

    def __getitem__(self,key):
        return self.data[key]

    def write(self, fout, key_re=None):
        for keyword in self.parsers:
            if (key_re is None or re.search(key_re, keyword)) and self.parsers[keyword][3] is not None:
                fout.write(self.parsers[keyword][3]())

    def write_atom_type_field(self):
        args_str = 'atom_type_field'
        args_str += ' {}'.format(self.data["atom_type_field"])
        return args_str+'\n'
    def parse_atom_type_field(self, args):
        args = self.parser_atom_type_field.parse_args(args)
        self.data["atom_type_field"] = args.field
        return "settings"

    def write_print_settings(self):
        return ''
    def parse_print_settings(self, args):
        args = self.parser_print_settings.parse_args(args)
        self.write(sys.stdout, key_re=args.keyword_regexp)
        return None

    def write_legend(self):
        args_str = 'legend '
        args_str += '-on' if self.data["legend"]['show'] else '-off'
        args_str += ' -position {} {}'.format(self.data["legend"]['position'][0],self.data["legend"]['position'][1])
        args_str += ' -spacing {}'.format(self.data["legend"]['spacing'])
        args_str += ' -sphere_scale {}'.format(self.data["legend"]['sphere_scale'])
        return args_str+'\n'
    def parse_legend(self, args):
        args = self.parser_legend.parse_args(args)

        got_setting = False
        if args.position is not None:
            self.data["legend"]['position'] = np.array(args.position)
            got_setting = True
        if args.offset is not None:
            self.data["legend"]['position'] += args.offset
            got_setting = True
        if args.spacing is not None:
            self.data["legend"]['spacing'] = args.spacing
            got_setting = True
        if args.sphere_scale is not None:
            self.data["legend"]['sphere_scale'] = args.sphere_scale
            got_setting = True


        if args.on:
            self.data["legend"]['show'] = True
        elif args.off:
            self.data["legend"]['show']= False
        elif not got_setting: # toggle
            self.data["legend"]['show'] = not self.data["legend"]['show']

        return "settings"

    def write_step(self):
        args_str = 'step'
        args_str += ' {}'.format(self.data["frame_step"])
        return args_str+'\n'
    def parse_step(self, args):
        args = self.parser_step.parse_args(args)
        self.data["frame_step"] = args.n
        return None

    def write_colormap (self):
        args_str = ''
        for colormap in self.data["colormaps"]:
            args_str += 'colormap {}'.format(colormap)
            for pt in self.data["colormaps"][colormap]:
                args_str += '   -P {} {} {} {}'.format(pt[0], pt[1], pt[2], pt[3])
            args_str += '\n'
        return args_str
    def parse_colormap(self, args):
        args = self.parser_colormap.parse_args(args)
        self.data["colormaps"][args.name] = args.colormap
        return "settings"

    def write_atom_type (self):
        args_str = ''
        all_data = self.data["atom_types"].get_all()
        for atom_type in all_data:
            args_str += 'atom_type {}'.format(atom_type)
            (color, colormap, radius, radius_field, material_args, bonding_radius) = all_data[atom_type]
            if color is not None:
                args_str += ' -color {} {} {}'.format(color[0], color[1], color[2])
            if colormap is not None:
                args_str += ' -colormap {} {}'.format(colormap[0], colormap[1])
            if radius is not None:
                args_str += ' -radius {}'.format(radius)
            if radius_field is not None:
                args_str += ' -radius_field {} {}'.format(radius_field[0], radius_field[1])
            if bonding_radius is not None:
                args_str += ' -bonding_radius {}'.format(bonding_radius)
            args_str += ' '+write_material_args(material_args)
            args_str += '\n'
        return args_str
    def parse_atom_type(self, args):
        args = self.parser_atom_type.parse_args(args)
        if args.radius_field is not None:
            args.radius_field[1] = float(args.radius_field[1])
        self.data["atom_types"].set_type(args.name, color=args.color, colormap=args.colormap, radius=args.radius, radius_field=args.radius_field,
                                         specular=args.specular, specular_radius=args.specular_radius, ambient=args.ambient, opacity=args.opacity,
                                         bonding_radius=args.bonding_radius)
        return "settings"

    def write_cell_box(self):
        args = 'cell_box {} {} {}'.format(self.data["cell_box"]["color"][0],
                                          self.data["cell_box"]["color"][1],
                                          self.data["cell_box"]["color"][2])
        args += ' -opacity {}'.format(self.data["cell_box"]["opacity"])
        args += ' -width {}'.format(self.data["cell_box"]["line_width"])
        return args+'\n'
    def parse_cell_box(self, args):
        args = self.parser_cell_box.parse_args(args)

        if args.color is not None:
            self.data["cell_box"]["color"] = (args.color[0], args.color[1], args.color[2])
            self.data["cell_box"]["prop"].SetColor(self.data["cell_box"]["color"])
        if args.opacity is not None:
            self.data["cell_box"]["opacity"] = args.opacity
            self.data["cell_box"]["prop"].SetOpacity(self.data["cell_box"]["opacity"])
        if args.width is not None:
            self.data["cell_box"]["line_width"] = args.width
            self.data["cell_box"]["prop"].SetLineWidth(self.data["cell_box"]["line_width"])
        return "color_only"

    def write_frame_label(self):
        args = 'frame_label'
        if self.data["frame_label"]["color"] is not None:
            args += ' -color {} {} {}'.format(self.data["frame_label"]["color"][0],
                                              self.data["frame_label"]["color"][1],
                                              self.data["frame_label"]["color"][2])
        if self.data["frame_label"]["fontsize"] is not None:
            args += ' -fontsize {}'.format(self.data["frame_label"]["fontsize"])
        if self.data["frame_label"]["string"] is not None:
            args += ' -string {}'.format(self.data["frame_label"]["string"])
        if self.data["frame_label"]["show"]:
            args += ' -on'
        else:
            args += ' -off'
        return args+'\n'
    def parse_frame_label(self, args):
        args = self.parser_frame_label.parse_args(args)

        got_setting = False
        if args.color is not None:
            self.data["frame_label"]["color"] = args.color
            got_setting = True
        if args.fontsize is not None:
            self.data["frame_label"]["fontsize"] = args.fontsize
            got_setting = True
        if args.string is not None:
            self.data["frame_label"]["string"] = " ".join(args.string)
            got_setting = True

        if args.on:
            self.data["frame_label"]["show"] = True
        elif args.off:
            self.data["frame_label"]["show"] = False
        elif not got_setting: # toggle
            self.data["frame_label"]["show"] = not self.data["frame_label"]["show"]
        self.data["frame_label"]["prop"].SetColor(self.data["frame_label"]["color"])
        self.data["frame_label"]["prop"].SetFontSize(self.data["frame_label"]["fontsize"])
        return "settings"

    def write_atom_label(self):
        args = 'atom_label'
        if self.data["atom_label"]["color"] is not None:
            args += ' -color {} {} {}'.format(self.data["atom_label"]["color"][0],
                                              self.data["atom_label"]["color"][1],
                                              self.data["atom_label"]["color"][2])
        if self.data["atom_label"]["fontsize"] is not None:
            args += ' -fontsize {}'.format(self.data["atom_label"]["fontsize"])
        if self.data["atom_label"]["string"] is not None:
            args += ' -string {}'.format(self.data["atom_label"]["string"])
        if self.data["atom_label"]["show"]:
            args += ' -on'
        else:
            args += ' -off'
        return args+'\n'
    def parse_atom_label(self, args):
        args = self.parser_atom_label.parse_args(args)

        got_setting = False
        if args.color is not None:
            self.data["atom_label"]["color"] = args.color
            got_setting = True
        if args.fontsize is not None:
            self.data["atom_label"]["fontsize"] = args.fontsize
            got_setting = True
        if args.string is not None:
            self.data["atom_label"]["string"] = args.string
            got_setting = True

        if args.on:
            self.data["atom_label"]["show"] = True
        elif args.off:
            self.data["atom_label"]["show"] = False
        elif not got_setting:
            self.data["atom_label"]["show"] = not self.data["atom_label"]["show"]

        self.data["atom_label"]["prop"].SetColor(self.data["atom_label"]["color"])
        self.data["atom_label"]["prop"].SetFontSize(self.data["atom_label"]["fontsize"])
        return "settings"

    def write_picked(self):
        args = 'picked {} {} {}'.format(self.data["picked"][0],
                                        self.data["picked"][1],
                                        self.data["picked"][2])
        return args+'\n'
    def parse_picked(self, args):
        args = self.parser_picked.parse_args(args)
        self.data["picked"]["color"] = (args.color[0], args.color[1], args.color[2])
        self.data["picked"]["prop"].SetColor(self.data["picked"]["color"])
        return "color_only"

    def write_background_color(self):
        args = 'background_color {} {} {}'.format(self.data["background_color"][0],
                                                  self.data["background_color"][1],
                                                  self.data["background_color"][2])
        return args+'\n'
    def parse_background_color(self, args):
        args = self.parser_background_color.parse_args(args)
        self.data["background_color"] = (args.color[0], args.color[1], args.color[2])
        return "color_only"
