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

class DavTKAtomTypes():
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
        sys.stderr.write("autogenerating {} with color {}\n".format(name, self.colors[self.autogen_used]))
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

        update_prop(self.types[name]["prop"], types.SimpleNamespace(**self.types[name]))

class DavTKSettings():
    parsers = {}

    _parser_atom_type_field = ThrowingArgumentParser(prog="atom_type_field", description="ASE at.arrays field to use for atom type")
    _parser_atom_type_field.add_argument("field", help="name of field ('Z' for atomic numbers, 'species' for chemical symbols)", default='Z')
    parsers["atom_type_field"] = _parser_atom_type_field

    _parser_print_settings = ThrowingArgumentParser(prog="print_settings", description="print settings")
    _parser_print_settings.add_argument("-keyword_regexp", type=str)
    parsers["print_settings"] = _parser_print_settings

    _parser_legend = ThrowingArgumentParser(prog="legend", description="control legend, toggle by default")
    _grp = _parser_legend.add_mutually_exclusive_group()
    _grp.add_argument("-on", action='store_true', help="enable legend")
    _grp.add_argument("-off", action='store_true', help="disable legend")
    _grp = _parser_legend.add_mutually_exclusive_group()
    _grp.add_argument("-position", type=int, nargs=2, help="position relative to bottom left corner of display"+
                                                         " (negative values relative to top right)", default=None)
    _grp.add_argument("-offset", type=int, nargs=2, help="offset relative to current position", default=None)
    _parser_legend.add_argument("-spacing", type=float, help="multiplier for spacing between rows", default=None)
    _parser_legend.add_argument("-sphere_scale", type=float, help="scaling factor for sphere radius", default=None)
    parsers["legend"] = _parser_legend

    _parser_step = ThrowingArgumentParser(prog="step", description="number of frames to skip in +/- and prev/next")
    _parser_step.add_argument("n", type=int, help="number of frames to step")
    parsers["step"] = _parser_step

    _parser_colormap = ThrowingArgumentParser(prog="colormap", description="repeated sequence of groups of 4 numbers: V R G B ...")
    _parser_colormap.add_argument("name", type=str)
    _parser_colormap.add_argument("-P", dest="colormap", nargs=4, action='append', type=float, metavar=('V', 'R', 'G', 'B'))
    parsers["colormap"] = _parser_colormap

    _parser_atom_type = ThrowingArgumentParser(prog="atom_type")
    _parser_atom_type.add_argument("name", type=str)
    _grp = _parser_atom_type.add_mutually_exclusive_group()
    _grp.add_argument("-color", "-c", nargs=3, type=float, default=None, metavar=("R", "G", "B"))
    _grp.add_argument("-colormap", nargs=2, default=None, metavar=("COLORMAP", "FIELD"))
    _grp = _parser_atom_type.add_mutually_exclusive_group()
    _grp.add_argument("-radius", "-rad", "-r", type=float, default=None)
    _grp.add_argument("-radius_field", nargs=2, metavar=("RADIUS_FIELD", "FACTOR"), default=None)
    _parser_atom_type.add_argument("-bonding_radius", type=float, default=None)
    add_material_args_to_parser(_parser_atom_type, "atom_type")
    parsers["atom_type"] = _parser_atom_type

    _parser_cell_box = ThrowingArgumentParser(prog="cell_box")
    _parser_cell_box.add_argument("-color", nargs=3, type=float, metavar=('R', 'G', 'B'), default=None)
    _parser_cell_box.add_argument("-opacity", type=float, default=None)
    _parser_cell_box.add_argument("-width", type=float, default=None)
    parsers["cell_box"] = _parser_cell_box

    _parser_picked = ThrowingArgumentParser(prog="picked")
    _parser_picked.add_argument("-color", nargs=3, type=float, metavar=('R', 'G', 'B'))
    parsers["picked"] = _parser_picked

    _parser_background_color = ThrowingArgumentParser(prog="background_color")
    _parser_background_color.add_argument("-color", nargs=3, type=float, metavar=('R', 'G', 'B'))
    parsers["background_color"] = _parser_background_color

    _parser_frame_label = ThrowingArgumentParser(prog="frame_label")
    _parser_frame_label.add_argument("-label", "-l", nargs='+',
        help="string, evaluating $( EXPRESSION ) and substituting ${STRING} with fields in atoms.info (or 'config_n'), or _NONE_", default=None)
    _parser_frame_label.add_argument("-color", "-c", nargs=3, type=float, default=None, metavar=("R", "G", "B"))
    _parser_frame_label.add_argument("-fontsize", type=int, default=None)
    _grp = _parser_frame_label.add_mutually_exclusive_group()
    _grp.add_argument("-on", action='store_true')
    _grp.add_argument("-off", action='store_true')
    parsers["frame_label"] = _parser_frame_label

    _parser_atom_label = ThrowingArgumentParser(prog="atom_label")
    _parser_atom_label.add_argument("-label", "-l", help="string to use for label, evaluating $( EXPRESSION ) and substituting $${STRING} with fields in atoms.arrays "+
                                                               "(or 'ID' for number of atom, 'Z' for atomic number, 'species' for chemical symbol), "+
                                                               "or '_NONE_'", default=None)
    _parser_atom_label.add_argument("-color", "-c", nargs=3, type=float, default=None, metavar=("R", "G", "B"))
    _parser_atom_label.add_argument("-fontsize", type=int, default=None)
    _grp = _parser_atom_label.add_mutually_exclusive_group()
    _grp.add_argument("-on", action='store_true')
    _grp.add_argument("-off", action='store_true')
    parsers["atom_label"] = _parser_atom_label

    def __init__(self):
        self.data = { "atom_types" : DavTKAtomTypes(), "colormaps" : {},
            "cell_box" : { "color" : [1.0, 1.0, 1.0], "opacity" : 1.0, "line_width" : 2.0, "prop" : None }, "background_color" : [0.0, 0.0, 0.0],
            "picked" : { "color" : [1.0, 1.0, 0.0], "opacity" : 1.0,  "prop" : None },
            "frame_label" : { "color" : [1.0, 1.0, 1.0], "fontsize" : 36, "prop" : None, "string" : "${config_n}", "show" : True },
            "atom_label" : { "color" : [1.0, 1.0, 1.0], "fontsize" : 24 , "prop" : None, "string" : "$${ID}", "show" : False},
            "frame_step" : 1, "legend" : { 'show' : False, 'position' : np.array([-10,-10]), 'spacing' : 1.0, 'sphere_scale' : 1.0 },
            'atom_type_field' : 'Z'
            }

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
            if (key_re is None or re.search(key_re, keyword)) and hasattr(self, "write_" + keyword):
                fout.write(getattr(self, "_write_" + keyword)())

    def _write_atom_type_field(self):
        args_str = 'atom_type_field'
        args_str += ' {}'.format(self.data["atom_type_field"])
        return args_str+'\n'
    def atom_type_field(self, field):
        self.data["atom_type_field"] = field
        return "settings"

    def _write_print_settings(self):
        return ''
    def print_settings(self, keyword_regexp):
        self.write(sys.stdout, key_re=keyword_regexp)
        return None

    def _write_legend(self):
        args_str = 'legend '
        args_str += '-on' if self.data["legend"]['show'] else '-off'
        args_str += ' -position {} {}'.format(self.data["legend"]['position'][0],self.data["legend"]['position'][1])
        args_str += ' -spacing {}'.format(self.data["legend"]['spacing'])
        args_str += ' -sphere_scale {}'.format(self.data["legend"]['sphere_scale'])
        return args_str+'\n'
    def legend(self, on=True, off=True, position=None, offset=None, spacing=None, sphere_scale=None):
        got_setting = False
        if position is not None:
            self.data["legend"]['position'] = np.array(position)
            got_setting = True
        if offset is not None:
            self.data["legend"]['position'] += offset
            got_setting = True
        if spacing is not None:
            self.data["legend"]['spacing'] = spacing
            got_setting = True
        if sphere_scale is not None:
            self.data["legend"]['sphere_scale'] = sphere_scale
            got_setting = True

        if on:
            self.data["legend"]['show'] = True
        elif off:
            self.data["legend"]['show']= False
        elif not got_setting: # toggle
            self.data["legend"]['show'] = not self.data["legend"]['show']

        return "settings"

    def _write_step(self):
        args_str = 'step'
        args_str += ' {}'.format(self.data["frame_step"])
        return args_str+'\n'
    def step(self, n):
        self.data["frame_step"] = n
        return None

    def _write_colormap(self):
        args_str = ''
        for colormap in self.data["colormaps"]:
            args_str += 'colormap {}'.format(colormap)
            for pt in self.data["colormaps"][colormap]:
                args_str += '   -P {} {} {} {}'.format(pt[0], pt[1], pt[2], pt[3])
            args_str += '\n'
        return args_str
    def colormap(self, name, colormap):
        self.data["colormaps"][name] = colormap
        return "settings"

    def _write_atom_type(self):
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
    def atom_type(self, name, color=None, colormap=None,
                      radius=None, radius_field=None, bonding_radius=None,
                      opacity=None, specular=None, specular_radius=None, ambient=None):
        if radius_field is not None:
            radius_field[1] = float(radius_field[1])
        self.data["atom_types"].set_type(name, color=color, colormap=colormap, radius=radius, radius_field=radius_field,
                                         specular=specular, specular_radius=specular_radius, ambient=ambient, opacity=opacity,
                                         bonding_radius=bonding_radius)
        return "settings"

    def _write_cell_box(self):
        args = 'cell_box -color {} {} {}'.format(self.data["cell_box"]["color"][0],
                                                 self.data["cell_box"]["color"][1],
                                                 self.data["cell_box"]["color"][2])
        args += ' -opacity {}'.format(self.data["cell_box"]["opacity"])
        args += ' -width {}'.format(self.data["cell_box"]["line_width"])
        return args+'\n'
    def cell_box(self, color=None, opacity=None, width=None):
        if color is not None:
            self.data["cell_box"]["color"] = (color[0], color[1], color[2])
            self.data["cell_box"]["prop"].SetColor(self.data["cell_box"]["color"])
        if opacity is not None:
            self.data["cell_box"]["opacity"] = opacity
            self.data["cell_box"]["prop"].SetOpacity(self.data["cell_box"]["opacity"])
        if width is not None:
            self.data["cell_box"]["line_width"] = width
            self.data["cell_box"]["prop"].SetLineWidth(self.data["cell_box"]["line_width"])
        return "color_only"

    def _write_picked(self):
        args = 'picked -color {} {} {}'.format(self.data["picked"]["color"][0],
                                               self.data["picked"]["color"][1],
                                               self.data["picked"]["color"][2])
        return args+'\n'
    def picked(self, color):
        self.data["picked"]["color"] = (color[0], color[1], color[2])
        self.data["picked"]["prop"].SetColor(self.data["picked"]["color"])
        return "color_only"

    def _write_background_color(self):
        args = 'background_color -color {} {} {}'.format(self.data["background_color"][0],
                                                         self.data["background_color"][1],
                                                         self.data["background_color"][2])
        return args+'\n'
    def background_color(self, color):
        self.data["background_color"] = (color[0], color[1], color[2])
        return "color_only"

    def _write_frame_label(self):
        args = 'frame_label'
        if self.data["frame_label"]["color"] is not None:
            args += ' -color {} {} {}'.format(self.data["frame_label"]["color"][0],
                                              self.data["frame_label"]["color"][1],
                                              self.data["frame_label"]["color"][2])
        if self.data["frame_label"]["fontsize"] is not None:
            args += ' -fontsize {}'.format(self.data["frame_label"]["fontsize"])
        if self.data["frame_label"]["string"] is not None:
            args += ' -label {}'.format(self.data["frame_label"]["string"])
        if self.data["frame_label"]["show"]:
            args += ' -on'
        else:
            args += ' -off'
        return args+'\n'
    def frame_label(self, label=None, color=None, fontsize=None, on=None, off=None):
        got_setting = False
        if color is not None:
            self.data["frame_label"]["color"] = color
            got_setting = True
        if fontsize is not None:
            self.data["frame_label"]["fontsize"] = fontsize
            got_setting = True
        if label is not None:
            self.data["frame_label"]["string"] = " ".join(label)
            got_setting = True

        if on:
            self.data["frame_label"]["show"] = True
        elif off:
            self.data["frame_label"]["show"] = False
        elif not got_setting: # toggle
            self.data["frame_label"]["show"] = not self.data["frame_label"]["show"]
        self.data["frame_label"]["prop"].SetColor(self.data["frame_label"]["color"])
        self.data["frame_label"]["prop"].SetFontSize(self.data["frame_label"]["fontsize"])
        return "settings"

    def _write_atom_label(self):
        args = 'atom_label'
        if self.data["atom_label"]["color"] is not None:
            args += ' -color {} {} {}'.format(self.data["atom_label"]["color"][0],
                                              self.data["atom_label"]["color"][1],
                                              self.data["atom_label"]["color"][2])
        if self.data["atom_label"]["fontsize"] is not None:
            args += ' -fontsize {}'.format(self.data["atom_label"]["fontsize"])
        if self.data["atom_label"]["string"] is not None:
            args += ' -label {}'.format(self.data["atom_label"]["string"])
        if self.data["atom_label"]["show"]:
            args += ' -on'
        else:
            args += ' -off'
        return args+'\n'
    def atom_label(self, label=None, color=None, fontsize=None, on=None, off=None):
        got_setting = False
        if color is not None:
            self.data["atom_label"]["color"] = color
            got_setting = True
        if fontsize is not None:
            self.data["atom_label"]["fontsize"] = fontsize
            got_setting = True
        if label is not None:
            self.data["atom_label"]["string"] = label
            got_setting = True

        if on:
            self.data["atom_label"]["show"] = True
        elif off:
            self.data["atom_label"]["show"] = False
        elif not got_setting:
            self.data["atom_label"]["show"] = not self.data["atom_label"]["show"]

        self.data["atom_label"]["prop"].SetColor(self.data["atom_label"]["color"])
        self.data["atom_label"]["prop"].SetFontSize(self.data["atom_label"]["fontsize"])
        return "settings"
