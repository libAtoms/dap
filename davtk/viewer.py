import os
import re, threading, queue, types
from argparse import ArgumentParser
import argparse

import warnings

from ase.atoms import Atoms
import ase.io

import vtk
from pathlib import Path
from davtk.settings import DavTKSettings
from davtk.state import *
from davtk.interactors import *

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt6.QtWidgets import QApplication, QMainWindow

from davtk.parse_utils import ThrowingArgumentParser, add_material_args_to_parser
from davtk.vtk_utils import new_prop, update_prop

try:
    import ffmpeg
    try:
        assert ffmpeg.input
    except AssertionError:
        warnings.warn("ffmpeg module loaded, but has no input function, maybe it's not actually ffmpeg-python")
        raise
except:
    ffmpeg = None

class Viewer(object):
    """Viewer of atomic configurations using VTK

    Parameters
    ----------
    at_list: list(ase.atoms.Atoms)
        list of atomic configurations
    win_size: (int, int), default (600, 600)
        initial window size
    win_name: str, default None
        name of window
    init_commands: list(str), default []
        list of commands to run before starting main loop (but after reading in settings and configs)
    """
    parsers = {}
    interactive = False
    _cmd_aliases = {}

    def __init__(self, at_list, win_size=(600, 600), win_name=None, init_commands=[]):
        if isinstance(at_list, Atoms):
            at_list = [at_list]
        try:
            for at in at_list:
                assert isinstance(at, Atoms)
                break
        except (AssertionError, TypeError) as exc:
            raise TypeError(f"at_list of type {type(at_list)} must be iterable of Atoms")

        # read settings from home dir and current dir
        self.davtk_state = None
        settings = DavTKSettings()

        if (Path.home() / ".daprc").exists():
            self._parse_file(Path.home() / ".daprc", settings=settings)
        if Path(".daprc").exists():
            self._parse_file(Path(".daprc"), settings=settings)

        if win_name is None:
            win_name = 'dap'

        # main window from Qt
        self.app = QApplication([win_name])
        self.renwin = QMainWindow()
        self.renwin.resize(*win_size)
        # interactor from Qt
        self.interactor = QVTKRenderWindowInteractor(self.renwin)
        self.renwin.setCentralWidget(self.interactor)
        self._application = False

        # A renderer
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(settings["background_color"])
        self.interactor.GetRenderWindow().AddRenderer(renderer)

        self.davtk_state = DaVTKState(at_list, settings, renderer, self.interactor)

        # add the custom styles for regular interaction and area selection
        sel_style = RubberbandSelect(self, parent=self.interactor)
        sel_style.SetDefaultRenderer(renderer)

        def_style = MouseInteractorHighLightActor(self, sel_style, parent=self.interactor)
        def_style.SetDefaultRenderer(renderer)
        # Switch to QTimer, as per
        #     https://stackoverflow.com/a/51608426
        # def_style.UseTimersOn()
        # self.interactor.CreateRepeatingTimer(100)

        self.interactor.SetInteractorStyle(def_style)

        # set up camera
        (min_pos, max_pos) = find_min_max(at_list)
        camera = renderer.GetActiveCamera()
        camera.ParallelProjectionOn()
        camera.SetParallelScale(np.max(max_pos-min_pos))
        extent = max(max_pos) - min(min_pos)
        camera.SetPosition([(max_pos[0]+min_pos[0])/2.0, -1000.0 + extent/2.0, (max_pos[2]+min_pos[2])/2.0])
        camera.SetViewUp(0.0, 0.0, 1.0)
        camera.SetFocalPoint((max_pos+min_pos)/2.0)
        camera.SetClippingRange(1000-2*extent, 1000+2*extent)

        l = vtk.vtkLight()
        l.SetPosition(100,70,200)
        l.SetLightTypeToCameraLight()
        renderer.GetLights().AddItem(l)

        # must do this rather late in the process
        ## renwin.SetWindowName(win_name)

        # read commands from atomic config headers
        for frame_i in range(len(self.davtk_state)):
            if "_vtk_commands" in at_list[frame_i].info:
                if not isinstance(at_list[frame_i].info["_vtk_commands"],str):
                    raise RuntimeError("Got _vtk_commands info field that is not a string")
                self.davtk_state.set_frame(str(frame_i))
                for cmd in at_list[frame_i].info["_vtk_commands"].split(";"):
                    try:
                        self._parse_line(cmd)
                    except Exception as exc:
                        print("Exception in _vtk_commands", exc)
                del at_list[frame_i].info["_vtk_commands"]

        self.davtk_state.set_frame("0")

        # now that atoms are read in and davtk_state exists, read any other commands (e.g. extra settings)
        for l in init_commands:
            for sub_l in l.split(";"):
                self._parse_line(sub_l)

        self.davtk_state.startup()

        self.davtk_state.prep_after_atoms_read()
        print("""DAP

        Use 'usage' for general usage info, and 'command -h' for detailed help on each command.
        Type 'h' in GUI window for GUI help message
        """)
        # Start
        self.renwin.show()
        self.interactor.Initialize()
        self.interactor.Start()

        self.davtk_state.update("cur")

        self.cmd_queue = queue.Queue(0)
        self.davtk_state.cmd_queue = self.cmd_queue


    def start(self, action, *args, **kwargs):
        """Start the viewer as the main application thread

        Parameters
        ----------
        action: function
            function to run in second thread, usually one that enqueues command line strings
        *args: list
            action positions args
        **kwargs: dict
            action keyword args
        """

        self._application = True
        self.cmd_thread = threading.Thread(target=action, args=args, kwargs=kwargs)
        self.cmd_thread.start()

        self.app.exec()


    def _parse_line(self, line, settings=None):
        """Read commands from a file and parse them

        Parameters
        ----------
        line: str
            line to parse
        settings: DavTKSettings, default None
            settings object (with parsers), default self.davtk_state.settings
        state: DavTKState / "_NONE_", default None
            davtk_state object (with parsers), default self.davtk_state, _NONE_ to not
            even try self.parsers which requires state (for parsing of settings before
            object is fully set up)
        """
        # settings and state from self if not explicitly specified
        if settings is None:
            settings = self.davtk_state.settings
        state = self.davtk_state

        if re.search(r'^\s*#', line) or len(line.strip()) == 0:
            # empty or comment
            return None

        args = line.strip().split()
        # special handling for command with special characters
        if re.search(r'[^a-zA-Z0-9_]', args[0]):
            # special character
            m = re.match(r'([^a-zA-Z0-9_]+)(.*)', args[0])
            if not m:
                raise ValueError(f'Command with special characters {args[0]} must start with them')
            args[0] = m.group(1)
            if len(m.group(2)) > 0:
                args.insert(1, m.group(2))

        cmd = self._cmd_aliases.get(args[0], args[0])

        matches_settings = [ k for k in settings.parsers.keys() if k.startswith(cmd) ]
        # only check for match in command parsers if state is actually available
        matches_cmds = [ k for k in self.parsers.keys() if k.startswith(cmd) ] if state is not None else []

        if len(matches_settings + matches_cmds) == 0:
            # no match
            raise ValueError("Unknown command '{}'\n".format(cmd))
        elif len(matches_settings + matches_cmds) > 1 and cmd not in matches_settings + matches_cmds:
            # more than 1 match and none are exact
            raise ValueError("Ambiguous command '{}', matches {}".format(cmd, matches_settings +
                                                                         matches_cmds))

        if cmd in matches_settings:
            # exact match in settings
            obj = settings
        elif cmd in matches_cmds:
            # exact match in cmds
            obj = self
        else:
            if len(matches_settings) == 1:
                obj = settings
                cmd = matches_settings[0]
            else:
                obj = self
                cmd = matches_cmds[0]

        parsed_args = obj.parsers[cmd].parse_args(args[1:])
        parsed_args = {k: v for k, v in vars(parsed_args).items() if v is not None}
        docstring = getattr(obj, cmd).__doc__
        if docstring is not None and "NEEDS_CMD_ARGS" in docstring:
            parsed_args["_cmd_args"] = args[1:]
        return getattr(obj, cmd)(**parsed_args)


    def _parse_file(self, filename, settings=None):
        """Read commands from a file and parse them

        Parameters
        ----------
        filename: str
            name of file to read from
        settings: DavTKSettings, default None
            settings object (with parsers), see parse_line()
        """
        refresh = {}

        with open(filename) as fin:
            for l in fin:
                refresh_status = self._parse_line(l, settings=settings)
                if refresh_status is not None:
                    refresh[refresh_status] = True

        if self.davtk_state is not None:
            for refresh_status in refresh:
                self.davtk_state.update(refresh_status)


    def cur_at(self):
        """Return currently viewed atoms object

        Returns
        -------
        cur_at: ase.atoms.Atoms cur atoms object
        """
        return self.davtk_state.cur_at()


    def _get_ats(self, all_frames):
        """Get a list of ats, either current or all

        Parameters
        ----------
        all_frames: bool
            return all atoms objects

        Returns
        -------
        ats: list(Atoms)
        """
        if all_frames:
            return self.davtk_state.at_list
        else:
            return [self.davtk_state.cur_at()]


    ####################################################################################################
    # commands and parsers
    ####################################################################################################

    _parser_usage = ThrowingArgumentParser(prog="usage", description="print usage message")
    _parser_usage.add_argument("-full_text", "-f", action="store_true", help="search full text of help message")
    _parser_usage.add_argument("glob", nargs="?", help="glob (shell style) to print usage for", default="*")
    parsers["usage"] = _parser_usage
    def usage(self, full_text=False, glob="*"):
        """Print usage message

        Parameters
        ----------
        full_text: bool, default False
            search for glob in full text of usage message
        glob: str, default "*"
            shell-style glob to match for keywords to print usage for

        Returns
        -------
        refresh: None
        """
        regexp = re.sub(r"\?", ".", glob)
        regexp = re.sub(r"\*", ".*", regexp)
        if not full_text:
            regexp = "^" + regexp + "$"

        for parsers, label in [(self.davtk_state.settings.parsers, "SETTINGS"),
                               (self.parsers, "COMMANDS")]:
            print(f"\n{label}:")
            for keyword in [k for k in sorted(parsers.keys()) if
                            (regexp is None or
                             re.search(regexp, k) or
                             (full_text and re.search(regexp, parsers[k].format_help(), re.MULTILINE)))]:
                print(keyword, parsers[keyword].format_usage(), end='')
        return None


    _parser_help = ThrowingArgumentParser(prog="help", description="print help message")
    _parser_help.add_argument("-full_text", "-f", action="store_true", help="search full text of help message")
    _parser_help.add_argument("glob", nargs="?", help="glob (shell style) to print usage for", default="*")
    parsers["help"] = _parser_help
    def help(self, full_text=False, glob="*"):
        """Print help message

        Parameters
        ----------
        full_text: bool, default False
            search for glob in full text of help message
        glob: str, default "*"
            shell-style glob to match for keywords to print usage for

        Returns
        -------
        refresh: None
        """
        regexp = re.sub(r"\?", ".", glob)
        regexp = re.sub(r"\*", ".*", regexp)
        if not full_text:
            regexp = "^" + regexp + "$"

        for parsers in (self.davtk_state.settings.parsers, self.parsers):
            for keyword in [k for k in sorted(parsers.keys()) if
                            (regexp is None or
                             re.match(regexp, k) or
                             (full_text and re.search(regexp, parsers[k].format_help(), re.MULTILINE)))]:
                print("--------------------------------------------------------------------------------")
                print(keyword, parsers[keyword].format_help(), end='')

        return None


    _parser_exit = ThrowingArgumentParser(prog="exit", description="exit viewer")
    parsers["exit"] = _parser_exit
    def exit(self):
        """Exit viewer

        Returns
        -------
        refresh: "exit"
        """
        if not self._application:
            for tl_widget in self.app.topLevelWidgets():
                tl_widget.close()
            for window in self.app.topLevelWindows():
                window.close()

        self.app.quit()
        return "exit"


    _parser_bang = ThrowingArgumentParser(prog="!", description="shell escape")
    _parser_bang.add_argument("command_line", nargs='+', help="shell escape command line")
    parsers["bang"] = _parser_bang
    _cmd_aliases["!"] = "bang"
    def bang(self, command_line):
        """Shell escape
        Aliases: !

        Returns
        -------
        refresh: None
        """
        os.system(" ".join(command_line))

        return None


    _parser_write_state = ThrowingArgumentParser(prog="write_state", description="write state of dap, including atomic configs, settings, and view")
    _parser_write_state.add_argument("-cur_frame_only", action="store_true")
    _parser_write_state.add_argument("filename", help="name of file to save into")
    parsers["write_state"] = _parser_write_state
    def write_state(self, filename, cur_frame_only=False):
        """write the state of the the viewer to a file

        Parameters
        ----------
        cur_frame_only: bool, default False
            only write the current frame, not all
        filename: str, required
            name of file to write to

        Returns
        -------
        refresh: None / "cur" / "settings" / "color_only" what things need to be refreshed
        """
        assert filename is not None, "filename is required"
        if not filename.endswith(".xyz") and not filename.endswith(".extxyz"):
            raise RuntimeError(f"Can only write state to extended xyz files, filename '{filename}' needs to end with .xyz or .extxyz")

        with open(filename+".settings", "w") as fout:
            self.davtk_state.settings.write(fout)
            fout.write("restore_view -data " + " ".join([str(v) for v in self.davtk_state.get_view()]) + "\n")

        with open(filename, "w") as fout:
            ats = [self.davtk_state.cur_at()] if cur_frame_only else self.davtk_state.at_list

            self.davtk_state.prep_for_atoms_write(ats)
            ase.io.write(fout, ats, format=ase.io.formats.filetype(filename, read=False))
        print("Settings written to '{0}.settings', read back with\n    dap -e \"read {0}.settings\" {0}".format(filename))
        return None


    _parser_restore_view = ThrowingArgumentParser(prog="restore_view",description="use stored viewing transform in ASE atoms object")
    _grp = _parser_restore_view.add_mutually_exclusive_group()
    _grp.add_argument("-name", help="name of info field to restore view from, default '_vtk_view' otherwise string will be prepended with '_vtk_view_'")
    _grp.add_argument("-data", type=float, nargs=12, help="list of 12 values defining view")
    parsers["restore_view"] = _parser_restore_view
    def restore_view(self, name=None, data=None):
        """Restore a view from a value saved to Atoms.info or raw numbers

        Parameters
        ----------
        name: str, default "_vtk_view"
            name of Atoms.info key to save view to, will be prepended by "_vtk_view_" if provided
        data: list(float)
            numbers specifying view, instead of array in Atoms.info

        Returns
        -------
        refresh: None
        """
        if data is None:
            # use name
            if name is None:
                name = "_vtk_view"
            else:
                name = "_vtk_view_" + name
        else:
            assert name is None, "At most one of name and data must be specified"

        if data:
            self.davtk_state.restore_view(data)
        else: # name
            at = self.davtk_state.cur_at()
            if name in at.info:
                self.davtk_state.restore_view(at.info[name])
            else:
                raise ValueError("restore_view info field '{}' not found".format(name))
        return None


    _parser_save_view = ThrowingArgumentParser(prog="save_view",description="store viewing transform in ASE atoms object")
    _parser_save_view.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _parser_save_view.add_argument("name", nargs="?", help="name of info field to save view into, default '_vtk_view', "
                                                           "otherwise string will be prepended with '_vtk_view_'")
    parsers["save_view"] = _parser_save_view
    def save_view(self, all_frames=False, name=None):
        """Save a view to an Atoms.info field

        Parameters
        ----------
        all_frames: bool, default False
            apply to all frames in list
        name: str, default "_vtk_view"
            name of Atoms.info key to save view to, will be prepended with "_vtk_view_" if specified
        """
        if name is None:
            name = "_vtk_view"
        else:
            name = "_vtk_view_" + name

        ats = self._get_ats(all_frames)

        view = self.davtk_state.get_view()
        for at in ats:
            at.info[name] = view

        return None


    _parser_movie = ThrowingArgumentParser(prog="movie",description="make a movie")
    _parser_movie.add_argument("-range", dest="frame_range", help="range of configs, in slice format start:[end+1]:[step]", default="::")
    _parser_movie.add_argument("-framerate",type=float,help="frames per second display rate", default=10.0)
    # _parser_movie.add_argument("-tmpdir",type=str,help="temporary directory for snapshots",default=".")
    _parser_movie.add_argument("-mag",type=int,help="magnification", default=1)
    _grp = _parser_movie.add_mutually_exclusive_group()
    _grp.add_argument("-ffmpeg_args",type=str,help="other ffmpeg args",default="-b:v 10M -pix_fmt yuv420p")
    _grp.add_argument("-raw_frames", dest="use_ffmpeg", action="store_false", help="Save raw frames without making into a movie with ffmpeg. "
                                                                                  "Requires that output_file contain a format-style substitution "
                                                                                  "for the frame number")
    _parser_movie.add_argument("output_file", help="Output file name, with ffmpeg-compatible suffix for movie (e.g. .mp4)"
                                                   "or frame-number format substitution (e.g. {:.03d}) and .png suffix for "
                                                   "individual frames")
    parsers["movie"] = _parser_movie
    def movie(self, output_file, frame_range="::", framerate=10.0, mag=1.0, ffmpeg_args="-b:v 10M -pix_fmt yuv420p",
                        use_ffmpeg=True):
        """Make a moview from frames

        Parameters
        ----------
        output_file: str
            filename to save movie to (requires ffmpeg), or format template with a single integer field for separate files for each frame
        frame_range: str, default "::"
            range of frames to use
        framerate: float, default 10
            frame rate in frames/sec
        mag: float, default 1.0
            magnification
        ffmpeg_args: str, default "-b:v 10M -pix_fmt yuv420p",
            arguments to pass ffmpeg
        use_ffmpeg: bool, default True
            whether to use ffmpeg or write to separate frames

        Returns
        -------
        refresh: None
        """
        m = re.search(r"^(\d*)(?::(\d*)(?::(\d*))?)?$", frame_range)
        if not m:
            raise SyntaxError("range '{}' is not in the expected format".format(frame_range))
        range_start = 0 if m.group(1) is None or len(m.group(1)) == 0 else int(m.group(1))
        range_end = len(self.davtk_state.at_list) if m.group(2) is None or len(m.group(2)) == 0 else int(m.group(2))
        range_interval = 1 if m.group(3) is None or len(m.group(3)) == 0 else int(m.group(3))
        frames = list(range(range_start, range_end, range_interval))

        if use_ffmpeg:
            if ffmpeg is None:
                raise RuntimeError('movie requires ffmpeg-python python module, unless -raw_frames is used')
        else:
            # check for format string in output_file
            try:
                if output_file.format(0) == output_file:
                    raise RuntimeError('No format string in -output_file "'+output_file+'"')
            except IndexError as exc:
                raise RuntimeError('Too many format strings in -output_file "'+output_file+'"') from exc

        for frame_i in frames:
            self.davtk_state.update(str(frame_i))

            if use_ffmpeg:
                data = self.davtk_state.snapshot(mag=mag).astype(np.uint8)
                data = np.ascontiguousarray(np.flip(data, axis=0))
                if frame_i == frames[0]:
                    # start ffmpeg
                    process = (
                        ffmpeg
                        .input('pipe:', format='rawvideo', pix_fmt='rgb24', s='{}x{}'.format(data.shape[1],data.shape[0]),
                               framerate=framerate)
                        .output(output_file, pix_fmt='yuv420p', framerate=framerate)
                        .overwrite_output()
                        .run_async(pipe_stdin=True)
                    )

                process.stdin.write(data)
            else:
                self.davtk_state.snapshot(output_file.format(frame_i), mag=mag)

        if use_ffmpeg:
            process.stdin.close()
            process.wait()

        # frames.append(frames[-1])
        # fmt_core = "0{}d".format(int(np.log10(len(frames)-1)+1))
        # py_fmt = ".{{:{}}}.png".format(fmt_core)
        # tmpfiles = []
        # with tempfile.NamedTemporaryFile(dir=tmpdir) as fout:
            # img_file_base = fout.name
            # for (img_i, frame_i) in enumerate(frames):
                # self.davtk_state.update(str(frame_i))
                # tmpfiles.append(img_file_base+py_fmt.format(img_i))
                # self.davtk_state.snapshot(tmpfiles[-1], 1)
        # print    ("ffmpeg -i {}.%{}.png -r {} {} {}".format(img_file_base, fmt_core, framerate, ffmpeg_args, "alt_"+output_file))
        # os.system("ffmpeg -i {}.%{}.png -r {} {} {}".format(img_file_base, fmt_core, framerate, ffmpeg_args, "alt_"+output_file))
        # for f in tmpfiles:
            # os.remove(f)

        return None


    _parser_frame = ThrowingArgumentParser(prog="frame", description="go to a particular frame")
    _parser_frame.add_argument("n", type=int, help="number of frame to go to")
    parsers["frame"] = _parser_frame
    _cmd_aliases["go"] = "frame"
    def frame(self, n):
        """Go to a specific frame
        Aliases: go

        Parameters
        ----------
        n: int / str
            number of frame to frame to, or offset if prefixed by "+" or "-"

        Returns
        -------
        refresh: None
        """
        self.davtk_state.update(str(n))

        return None


    _parser_next = ThrowingArgumentParser(prog="next", description="go forward a number of frames")
    _parser_next.add_argument("n", type=int, nargs='?', help="number of frames to increment (default set by 'step' command)")
    parsers["next"] = _parser_next
    _cmd_aliases["n"] = "next"
    _cmd_aliases["+"] = "next"
    def next(self, n=None):
        """Go to a next (one or more) frame
        Aliases: n, +

        Parameters
        ----------
        n: int, default None
            number of frame to increment by, default to number set by 'step' command

        Returns
        -------
        refresh: None
        """
        assert n is None or n > 0, "n must be > 0"
        self.davtk_state.update("+{}".format(n if n is not None else self.davtk_state.settings["frame_step"]))
        return None


    _parser_previous = ThrowingArgumentParser(prog="previous", description="go forward a number of frames")
    _parser_previous.add_argument("n", type=int, nargs='?', help="number of frames to decrement (default set by 'step' command)")
    parsers["previous"] = _parser_previous
    _cmd_aliases["p"] = "previous"
    _cmd_aliases["-"] = "previous"
    def previous(self, n=None):
        """Go to a previous (one or more) frame
        Aliases: p, -

        Parameters
        ----------
        n: int, default None
            number of frame to decerement by, default to number set by 'step' command

        Returns
        -------
        refresh: None
        """
        assert n is None or n > 0, "n must be > 0"
        self.davtk_state.update("-{}".format(n if n is not None else self.davtk_state.settings["frame_step"]))
        return None


    _parser_pick = ThrowingArgumentParser(prog="pick", description="pick atom(s) by ID")
    _parser_pick.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _parser_pick.add_argument("n", type=int, nargs='+')
    parsers["pick"] = _parser_pick
    def pick(self, n, all_frames=False):
        """Pick atoms by their number

        Parameters
        ----------
        n: list(int)
            indices of atoms to pick
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        ats = self._get_ats(all_frames)

        for at in ats:
            if "_vtk_picked" not in at.arrays:
                at.new_array("_vtk_picked",np.array([False]*len(at)))
            at.arrays["_vtk_picked"][n] = True
        return "cur"


    _parser_unpick = ThrowingArgumentParser(prog="unpick", description="unpick picked atoms and/or bonds (default both)")
    _parser_unpick.add_argument("-all_frames", action="store_true",help="apply to all frames")
    _parser_unpick.add_argument("-atoms", action="store_true")
    _parser_unpick.add_argument("-bonds", action="store_true")
    parsers["unpick"] = _parser_unpick
    def unpick(self, atoms=False, bonds=False, all_frames=False):
        """Unick atoms and/or bonds

        Unpicks all if atoms and bonds are both False

        Parameters
        ----------
        atoms: bool, default False
            unpick atoms (only atoms if True and bonds is False)
        bonds: bool, default False
            unpick bonds (only bonds if True and atoms is False)
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        ats = self._get_ats(all_frames)

        none_selected = not any([atoms, bonds])
        for at in ats:
            if none_selected or atoms:
                if "_vtk_picked" not in at.arrays:
                    continue
                at.arrays["_vtk_picked"][:] = False
            if none_selected or bonds:
                if hasattr(at, "bonds"):
                    for i_at in range(len(at)):
                        for b in at.bonds[i_at]:
                            b["picked"] = False
        return "cur"


    _parser_delete = ThrowingArgumentParser(prog="delete", description="delete objects (picked by default)")
    _parser_delete.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _parser_delete.add_argument("-atom_index", type=int, nargs='+', help="delete atoms by these indices ", default=None)
    _parser_delete.add_argument("-bond_name", nargs='+', help="delete bonds with these names ", default=None)
    parsers["delete"] = _parser_delete
    def delete(self, atom_index=None, bond_name=None, all_frames=False):
        """delete atoms and/or bonds (only picked by default)

        Parameters
        ----------
        atoms_index: list(int), default None
            list of indices of atoms to delete
        bond_name: list(str)
            list of bond names to delete
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        ats = self._get_ats(all_frames)

        for at in ats:
            did_something = False
            if atom_index is not None:
                self.davtk_state.delete(at, atom_selection=atom_index)
                did_something = True
            if bond_name is not None:
                for name in bond_name:
                    self.davtk_state.delete(at, bond_selection=name)
                did_something = True
            if not did_something:
                self.davtk_state.delete(at, atom_selection="picked", bond_selection="picked")

        return "cur"

    _parser_supercell = ThrowingArgumentParser(prog="supercell",description="create supercell of current cell")
    _parser_supercell.add_argument("-all_frames",action="store_true",help="apply to all frames")
    _parser_supercell.add_argument("-wrap",type=bool,default=None,help="wrap position into cell after creating supercell")
    _parser_supercell.add_argument("n", type=int, nargs='+',
                                   help="number of periodic images of cell to create: scalar (factor in all 3 dirs) or "
                                        "3-vector (factor in each dir) or 9 ints with 3 x 3 transformation matrix of new cell")
    parsers["supercell"] = _parser_supercell
    def supercell(self, n, wrap=None, all_frames=False):
        """Create supercell of configuration(s)

        Parameters
        ----------
        n: int / list(int)
            scalar, vector, or 3x3 matrix describing supercell
        wrap: bool, default False
            wrap atoms after creating supercell, default False for scalar or vector n and True for 3x3 matrix n
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        if all_frames:
            frame_list=None
        else:
            frame_list="cur"

        # promote to list
        try:
            _ = len(n)
        except TypeError:
            n = [n]

        assert len(n) in (1, 3, 9), f"Allowed lengths for n 1, 3, 9, got n = {n}"

        if wrap is None:
            wrap = False

        if len(n) == 1:
            n = n * 3
        elif len(n) == 9:
            if wrap is None:
                wrap = True

        self.davtk_state.supercell(n, wrap, frames=frame_list)
        return None


    _parser_images = ThrowingArgumentParser(prog="images", description="show images of cell")
    _parser_images.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _parser_images.add_argument("-lattice_field", help="Field to use for lattice, _CART_ for Cartesian axes", default='_CELL_')
    _parser_images.add_argument("r", type=float, nargs='+', help="range in lattice coordinates (floating point) of images to display. "
        "If scalar or 3-vector, padding in addition to original cell (i.e. -R -- 1+R).  If 6-vector, r0_min -- r0_max, "
        "r1_min -- r1_max, r2_min -- r2_max.  If scalar < 0, reset")
    parsers["images"] = _parser_images
    def images(self, r, lattice_field="_CELL_", all_frames=False):
        """Display periodic images without varying cell box

        Parameters
        ----------
        r: float / list(float) / None
            if scalar or len 1 or 3, range beyond original (i.e. -r .. 1+r). If len 6, min and max
            in each direction. If None, reset
        lattice_field: str / "_CELL_" / "_CART_" default "_CELL"
            vectors to use when generating extra images, _CELL_ for Atoms.cell, _CART_ for Cartesian
            unit vectors, or string for cell stored in Atoms.info key
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        ats = self._get_ats(all_frames)

        if r is not None:
            try:
                _ = len(r)
            except TypeError:
                r = [r]

        assert len(r) in (1, 3, 6), f"r must be scalar, 3-vector, or 3 min, max pairs, got {r}"

        if r is not None:
            if len(r) == 1:
                if r[0] < 0:
                    r = None
                else:
                    r = [-r[0], 1.0 + r[0]] * 3
            elif len(r) == 3:
                r = [-r[0], 1.0 + r[0], -r[1], 1.0 + r[1], -r[2], 1.0 + r[2]]

        for at in ats:
            if r is None:
                at.info.pop("_vtk_images", None)
                at.info.pop("_vtk_images_cell_field", None)
            else:
                at.info["_vtk_images"] = r
                at.info["_vtk_images_cell_field"] = lattice_field

        return None


    _parser_vectors = ThrowingArgumentParser(prog="vectors", description="Draw vectors")
    _parser_vectors.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _parser_vectors.add_argument("-field", type=str, help="atom field to use for vectors (scalar or 3-vector)", default=None)
    _parser_vectors.add_argument("-color", type=str, metavar=["COLOR_SCHEME"], nargs='+',
                                 help="The string 'atom' (color by atom), or one color: R G B, "
                                      "or two colors : RUP GUP BUP RDOWN GDOWN BDOWN", default=None)
    _parser_vectors.add_argument("-radius", type=float, help="cylinder radius", default=None)
    _parser_vectors.add_argument("-scale", type=float, help="scaling factor from field value to cylinder length", default=None)
    _parser_vectors.add_argument("-delete", action='store_true', help="disable")
    parsers["vectors"] = _parser_vectors
    def vectors(self, field=None, color=None, radius=None, scale=None, delete=False, all_frames=False):
        """Display vectors from each atom

        Parameters
        ----------
        field: str, default "magmoms"
            Atoms.arrays key (either scalar of 3-vector) or special keys "forces" or "magmoms" or "initial_magmoms"
        color: "atom" / (v, v, v) / (v, v, v, v, v, v) v is str or float
            color, either by atom or a color triplet or two triplets one for up and one for down (for scalar "vector" properties)
        radius: float
            radius of vector cylinder
        scale: float
            scale between vectors in raw units and positions units
        delete: bool
            delete bonds in field
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        if delete:
            assert field is None and color is None and radius is None and scale is None, "Got delete, no other arguments may be specified"

        if color is not None and isinstance(color, str):
            color = [color]

        sign_colors = None
        if color is not None:
            assert len(color) in (1, 3, 6), f"Color not length 1, 3, or 6, got {color}"
            if len(color) == 1: # color by atom
                if color[0] != "atom":
                    raise ValueError("Got single word color argument '{}', not 'atom'".format(color))
                color = "atom"
            elif len(color) == 3: # single color
                try:
                    color = [float(v) for v in color]
                except:
                    raise ValueError("Got 3-word color '{}' not convertible to float".format(color))
            else: # color by sign
                try:
                    sign_colors = [float(v) for v in color]
                except:
                    raise ValueError("Got 3-word color '{}' not convertible to float".format(color))
                color = "sign"

        ats = self._get_ats(all_frames)

        for at in ats:
            if delete:
                del at.info["_vtk_vectors"]
            else:
                if "_vtk_vectors" not in at.info:
                    at.info["_vtk_vectors"] = { "field" : "magmoms", "color" : "atom", "radius" : 0.1, "scale" : 1.0,
                                                "sign_colors" : [1.0, 0.0, 0.0,    0.0, 0.0, 1.0] }
                for p in ["field", "color", "sign_colors", "radius", "scale" ]:
                    if locals()[p] is not None:
                        at.info["_vtk_vectors"][p] = locals()[p]
        return None

    _parser_bond = ThrowingArgumentParser(prog="bond",description="Create bonds")
    _parser_bond.add_argument("-all_frames",action="store_true",help="apply to all frames")
    _parser_bond.add_argument("-name",type=str,help="name of bond set", default="default")
    _parser_bond.add_argument("-T",type=str,help="Restrict one member to given atom_type", default="*")
    _parser_bond.add_argument("-Tn",type=str,help="Restrict other member to given atom_type", default="*")
    _grp = _parser_bond.add_mutually_exclusive_group()
    _grp.add_argument("-picked", action='store_true', help="bond a pair of picked atoms with MIC")
    _grp.add_argument("-n", action='store', type=int, nargs=2, help="bond a pair of atoms with given IDs", default=None)
    _grp.add_argument("-cutoff", "-cut", "-c", action='store', type=float, nargs='+', help="bond atoms with max (one argument) or min-max (two arguments) cutoffs", default=None)
    _grp.add_argument("-auto", action='store_true', help="bond by previously set per-atom bonding radius")
    _parser_bond.add_argument("-radius", type=float, help="radius of bond cylinders", default=None)
    _parser_bond.add_argument("-color", nargs=3, type=float, metavar=("R","G","B"), help="color for bonds", default=None)
    add_material_args_to_parser(_parser_bond, "bonds")
    _grp = _parser_bond.add_mutually_exclusive_group()
    _grp.add_argument("-delete", action='store_true', help="delete existing bond")
    _grp.add_argument("-list_bonds", action='store_true', help="list existing bond")
    _parser_bond.add_argument("-no_pbc", action='store_true', help="no bonds across periodic boundary conditions")
    parsers["bond"] = _parser_bond
    def bond(self, name="default", T="*", Tn="*", picked=False, n=None, cutoff=None, auto=False,
                 radius=None, color=None, opacity=None, specular=None, specular_radius=None, ambient=None,
                 delete=False, list_bonds=False, no_pbc=False, all_frames=False):
        """Create bonds between atoms

        Parameters
        ----------
        name: str, default "default"
            name to attach to set of bonds
        T: str, default "*"
            atom type for one end of bond
        Tn: str, default "*"
            atom type for other end of bond
        picked: bool
            bond all atoms that are picked
        n: list(int)
            ids of atoms to bond
        cutoff: float
            bond all atoms within cutoff
        auto: bool
            bond all atoms within per-atom-type radius
        radius: float
            radius of bond cylinder
        color: (float, float, float)
            color of bond
        opacity: float
            opacity of bond
        specular: float
            amount of specular reflection
        specular_radius: float
            shininess of specular reflection
        ambient: float
            amount of ambient color
        delete: bool
            delete bonds with specified name
        list_bonds: bool
            list bond sets that exist
        no_pbc: bool
            do not draw bonds across pbcs
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        ats = self._get_ats(all_frames)
        davtk_state = self.davtk_state

        creating_bonds = picked or (n is not None) or (cutoff is not None) or auto

        if creating_bonds:
            # check for conflicting args
            if delete or list_bonds:
                raise RuntimeError("can't create bonds and also -delete or -list_bonds")
            if cutoff is not None and len(cutoff) > 2:
                raise ValueError("bond got -cutoff with {} > 2 values".format(len(cutoff)))

        # create new property if needed (even if not creating bonds)
        if name not in davtk_state.bond_prop:
            prop = new_prop( types.SimpleNamespace( color = (1.0, 1.0, 1.0), opacity = 1.0,
                specular = 0.0, specular_radius = 0.1, ambient = 0.2 ) )
            prop.radius = 0.2
            davtk_state.bond_prop[name] = prop

        if not delete and not list_bonds:
            # update property
            update_prop(davtk_state.bond_prop[name], types.SimpleNamespace(color=color, opacity=opacity,
                        specular=specular, specular_radius=specular_radius, ambient=ambient))
            if radius is not None:
                davtk_state.bond_prop[name].radius = radius

        if creating_bonds:
            # actually create bonds
            for at in ats:
                if picked:
                    if T != '*' or Tn != '*':
                        raise RuntimeError("bond by picked can't specify -T or -Tn")
                    davtk_state.bond(at, name, T, Tn, "picked")
                elif n:
                    if T != '*' or Tn != '*':
                        raise RuntimeError("bond by n (indices) can't specify -T or -Tn")
                    davtk_state.bond(at, name, T, Tn, ("n", n))
                elif cutoff:
                    davtk_state.bond(at, name, T, Tn, ("cutoff", cutoff), no_pbc=no_pbc)
                elif auto:
                    davtk_state.bond(at, name, T, Tn, ("cutoff", None), no_pbc=no_pbc)
        else: # not creating, either delete or list or just modifying an existing name
            if delete:
                for at in ats:
                    davtk_state.delete(at, atom_selection=None, bond_selection=name)
            elif list_bonds:
                if len(davtk_state.bond_prop) > 0:
                    print("defined (not necessarily used) bond names: ", list(davtk_state.bond_prop.keys()))
                else:
                    print("no bond names defined")
            # else: just modifying existing proprerty

        # TODO
        # look for unused bond names and remove their props

        if creating_bonds or delete:
            return "cur"
        elif radius is None:
            return "color_only"
        else:
            return "settings"


    _parser_snapshot = ThrowingArgumentParser(prog="snapshot",description="write snapshot")
    _parser_snapshot.add_argument("-frame_range", help="slice (start:end:step) to generate snapshots for", default=None)
    _parser_snapshot.add_argument("-mag", type=int, help="magnification", default=1)
    _parser_snapshot.add_argument("file", type=str, help="filename (png format), with single format() {} string for "
                                                         "int argument if -frame_range is not None")
    parsers["snapshot"] = _parser_snapshot
    def snapshot(self, file, frame_range=None, mag=1.0):
        """Take a snapshot of one or more frames

        Parameters
        ----------
        file: str
            name of file to save to, must have single format string for int frame number of frame_range is not None
        mag: float, default 1
            magnification relative to current window size
        frame_range: str
            range (start:end:step) to take snapshots of, default cur frame only
        """
        if frame_range is not None:
            if file.format(1) == file:
                raise RuntimeError('Need format string in file argument when slice is not None')
            range_slice = [int(s) if len(s) > 0 else None for s in frame_range.split(':')]
            if len(range_slice) > 3:
                raise RuntimeError('-slice format not [start[:[end][:step]]]')
            range_slice += [None] * (3-len(range_slice))
            for frame_i in range(len(self.davtk_state.at_list))[slice(*range_slice)]:
                self.davtk_state.update(str(frame_i))
                self.davtk_state.snapshot(file.format(frame_i), mag)
        else:
            self.davtk_state.snapshot(file, mag)

        return None


    _parser_X = ThrowingArgumentParser(prog="X",description="execute python code (Atoms object available as 'atoms', "
        "DavTKSettings as 'settings'). Use '--' to separate command from X's arguments")
    _parser_X.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _parser_X.add_argument("-global_context", "-g", action="store_true", help="run in global scope (for 'import', e.g.)")
    _parser_X.add_argument("command_line", nargs='+', help="X command line words")
    parsers["X"] = _parser_X
    def X(self, command_line, all_frames=False, global_context=False):
        ase_command = " ".join(command_line)
        ats = self._get_ats(all_frames)
        settings = self.davtk_state.settings

        for atoms in ats:
            globals()["atoms"] = atoms
            if global_context:
                exec(ase_command) in globals(), globals()
            else:
                exec(ase_command)

        return "cur"


    _parser_read = ThrowingArgumentParser(prog="read",description="read commands from file(s)")
    _parser_read.add_argument("filename", nargs='+', help="filenames")
    parsers["read"] = _parser_read
    def read(self, filename):
        """Read commands from a file

        Parameters
        ----------
        filename: str / list(str)
            one or more files to read from
        """
        if isinstance(filename, str):
            filename = [filename]

        for f in filename:
            self._parse_file(f)


    _parser_override_frame_label = ThrowingArgumentParser(prog="override_frame_label", description="set per-frame label")
    _parser_override_frame_label.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _parser_override_frame_label.add_argument("label", help="per-configuration string, substituting ${FIELD} with fields "+
                                                            "in atoms.info (or 'config_n'), or _NONE_")
    parsers["override_frame_label"] = _parser_override_frame_label
    def override_frame_label(self, label=None, all_frames=False, func=None):
        """Set non-default frame label

        Parameters
        ----------
        label: str, default None
            label string, substituting ${FIELD} with atoms.info or 'config_n' for config_n, or "_NONE_" for no label
        all_frames: bool, default False
            apply to all frames rather than just current
        func: callable, default None
            callable to create label, must take two arguments, Atoms and int at_i, and return a string
        """

        ats = self._get_ats(all_frames)

        assert sum([label is None, func is None]) == 1, "Must get exactly one of label and func"

        for at_i, at in enumerate(ats):
            if label is not None:
                at.info["_vtk_frame_label_string"] = label
            else:
                at.info["_vtk_frame_label_string"] = func(at, at_i)

        return "cur"

    _parser_override_atom_label = ThrowingArgumentParser(prog="override_atom_label", description="override atom label on a per-atom basis")
    _parser_override_atom_label.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _grp = _parser_override_atom_label.add_mutually_exclusive_group()
    _grp.add_argument("-picked", action='store_true', help="label picked atoms with STRING, default", default=None)
    _grp.add_argument("-n", help="comma separate lst of atom ids to label", default=None)
    _grp.add_argument("-mask", nargs='+', help="expression that evaluates "+
                                               "to list of indices or logical array with dim len(atoms)", default=None)
    _parser_override_atom_label.add_argument("label", nargs='+', help="string to set as per-atom label", default=[])
    parsers["override_atom_label"] = _parser_override_atom_label
    def override_atom_label(self, picked=None, n=None, mask=None, label=None, all_frames=False, func=None):
        """override atom label

        Parameters
        ----------
        picked: bool
            apply to picked atoms
        n: int / list(int) / str
            apply to list of atom ids, str is interpreted as comma separated list
        mask: str
            apply to atoms selected by mask from evaluating mask expression (list of indices or list of bools)
        label: str
            content of label
        all_frames: bool, default False
            apply to all frames rather than just current
        func: callable, default None
            callable to create label, must single Atoms argument return an array of strings, one for each atom
        """
        ats = self._get_ats(all_frames)

        for at in ats:
            if func is not None:
                assert label is None, "Can only accept one of label or func"
                new_label = np.asarray(func(at))
            else:
                new_label = np.asarray([" ".join(label)] * len(at))
            if "_vtk_label" not in at.arrays or len(new_label) > len(at.arrays["_vtk_label"][0]):
                a = np.asarray([" " * max([len(l) for l in new_label])] * len(at))
                if "_vtk_label" in at.arrays:
                    a[:] = at.arrays["_vtk_label"]
                    del at.arrays["_vtk_label"]
                at.new_array("_vtk_label", a)
            if n:
                for range_item in n.split(","):
                    exec(f'at.arrays["_vtk_label"][{range_item}] = new_label[{range_item}]')
            elif mask:
                atoms = at
                expr = " ".join(mask)
                exec(f'atoms.arrays["_vtk_label"][{expr}] = new_label[{expr}]')
            else: # -picked, or nothing selected (default)
                if "_vtk_picked" not in at.arrays:
                    return None
                for i_at in range(len(at)):
                    if at.arrays["_vtk_picked"][i_at]:
                        at.arrays["_vtk_label"][i_at] = new_label[i_at]
        return "cur"


    _parser_measure = ThrowingArgumentParser(prog="measure",description="measure some quantities for picked objects (default) or listed atoms")
    _parser_measure.add_argument("-all_frames",action="store_true",help="apply to all frames")
    _parser_measure.add_argument("-n",action="store",type=int,nargs='+',help="ID numbers of atoms to measure", default = None)
    parsers["measure"] = _parser_measure
    def measure(self, n=None, all_frames=False):
        """measure positions/distances/angles of atoms and print to stdout

        Parameters
        ----------
        n: int / list(int)
            indices of atoms to measure, default picked
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        if n is not None:
            try:
                _ = len(n)
            except TypeError:
                n = [n]

        if all_frames:
            frames = range(len(self.davtk_state.at_list))
        else:
            frames = [self.davtk_state.cur_frame]

        for frame_i in frames:
            print("Frame:",frame_i)
            self.davtk_state.measure(n, frame_i)
        return None


    def _polyhedra(self, name, T, Tn, color, opacity, specular, specular_radius, ambient,
                   cutoff, bond_name, indices, delete, list_polys, all_frames):

        assert not (indices is not None and (T is not None or Tn is not None or cutoff is not None or bond_name is not None))

        ats = self._get_ats(all_frames)
        davtk_state = self.davtk_state

        creating_polyhedra = cutoff is not None or bond_name is not None or indices is not None
        modifying_polyhedra = not (creating_polyhedra or delete or list_polys)

        polyhedra_names = list(set([n.replace("_vtk_polyhedra_","") for at in ats for n in at.arrays if n.startswith("_vtk_polyhedra")]))

        # create new property if needed
        if creating_polyhedra and name not in davtk_state.polyhedra_prop:
            prop = new_prop( types.SimpleNamespace( color = (0.5, 0.5, 1.0), opacity = 0.5,
                specular = 0.7, specular_radius = 0.1, ambient = 0.1 ) )
            davtk_state.polyhedra_prop[name] = prop
        # update to desired values
        if creating_polyhedra or modifying_polyhedra:
            if name not in davtk_state.polyhedra_prop:
                sys.stderr.write(f"ERROR: can't modify {name} which is not a known polyhedra name {polyhedra_names}")
            # update property
            update_prop(davtk_state.polyhedra_prop[name], types.SimpleNamespace(color=color, opacity=opacity,
                        specular=specular, specular_radius=specular_radius, ambient=ambient))

        if creating_polyhedra:
            # actualy create
            if indices is not None:
                for at in ats:
                    davtk_state.arb_polyhedra(at, name, indices)
            else:
                # check for required args (cutoff/bond_name already checked above)
                if T is None:
                    raise RuntimeError("polyhedra requires -T to create new polyhedra")
                for at in ats:
                    davtk_state.coordination_polyhedra(at, name, T, Tn, cutoff, bond_name)
        elif delete is not None:
            for at in ats:
                if "_vtk_polyhedra_" + delete in at.arrays:
                    del at.arrays["_vtk_polyhedra_" + delete]
            # TODO: look for unused polyhedra names and remove their props ?
        elif list_polys:
            print("polyhedra sets:", polyhedra_names)
            return
        # else: just modifying existing property

        if creating_polyhedra or delete:
            return "cur"
        else:
            return "color_only"


    # Logic of which arguments can coexist is way too messy here. Not sure how to fix.
    _parser_polyhedra = ThrowingArgumentParser(prog="polyhedra",description="draw coordination polyhedra")
    _parser_polyhedra.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _parser_polyhedra.add_argument("-name", help="name of polyhedron set", default="default")
    _parser_polyhedra.add_argument("-T", type=str, help="atom_type for polyhedron center, required for creating",  default=None)
    _parser_polyhedra.add_argument("-Tn", type=str, help="atom_type for polyhedron neighbors",  default=None)
    _parser_polyhedra.add_argument("-color", nargs=3, type=float, metavar=("R","G","B"), help="color for polyhedra", default=None)
    add_material_args_to_parser(_parser_polyhedra, "polyhedra")
    _grp = _parser_polyhedra.add_mutually_exclusive_group()
    _grp.add_argument("-cutoff", type=float, help="center-neighbor cutoff distance",  default=None)
    _grp.add_argument("-bond_name", help="name of existing bond to use for neighbors", default=None)
    _grp.add_argument("-delete", action='store', metavar="NAME", help="name of polyhedra to delete")
    _grp.add_argument("-list_polys", action='store_true', help="list existing polyhedra")
    parsers["polyhedra"] = _parser_polyhedra
    def polyhedra(self, name=None, T=None, Tn=None, color=None, opacity=None, specular=None, specular_radius=None, ambient=None,
                      cutoff=None, bond_name=None, delete=None, list_polys=None, all_frames=False):
        """Create coordination polyhedra based on bonds

        Parameters
        ----------
        name: str
            name to associate with polyhedra
        T: str
            atom type of center of polyhedra
        Tn: str
            atom type of corners of polyhedra, defaults to any atom
        color: (float, float, float)
            polyhedron color
        opacity: float
            polyhedron opacity
        specular: float
            polyhedron specular reflection
        specular_radius: float
            polyhedron specular shininess
        ambient: float
            polyhedron ambient light
        cutoff: float
            define polyhedron based on neighbors within this cutoff
        bond_name: str
            define polyhedron based on bonds with this name
        delete: str
            delete polyhedra with this name
        list_polys: bool
            list existing polyhedra
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        return self._polyhedra(name, T, Tn, color, opacity, specular, specular_radius, ambient,
                               cutoff, bond_name, None, delete, list_polys, all_frames)


    _parser_arb_polyhedra = ThrowingArgumentParser(prog="arb_polyhedra",description="draw arbitrary polyhedra connecting listed atoms")
    _parser_arb_polyhedra.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _parser_arb_polyhedra.add_argument("-name", help="name of polyhedron set, for later reference with -delete", default="default")
    _parser_arb_polyhedra.add_argument("-color", nargs=3, type=float, metavar=("R","G","B"), help="color for polyhedra", default=None)
    add_material_args_to_parser(_parser_arb_polyhedra, "polyhedra")
    _grp = _parser_arb_polyhedra.add_mutually_exclusive_group()
    _grp.add_argument("-indices", type=int, action='append', nargs='+', help="atomic indices of polyhedron", default=None)
    _grp.add_argument("-delete", action='store', metavar="NAME", help="name of polyhedra to delete")
    _grp.add_argument("-list_polys", action='store_true', help="list existing polyhedra")
    parsers["arb_polyhedra"] = _parser_arb_polyhedra
    def arb_polyhedra(self, name=None, color=None, opacity=None, specular=None, specular_radius=None, ambient=None,
                          indices=None, delete=None, list_polys=None, all_frames=False):
        """Create arbitrary polyhedra

        Parameters
        ----------
        name: str
            name associated with these polyhedra
        color: (float, float, float)
            color of polyhedra
        opacity: float
            polyhedron opacity
        specular: float
            polyhedron specular reflection
        specular_radius: float
            polyhedron specular shininess
        ambient: float
            polyhedron ambient light
        indices: list(int)
            indices of atoms defining polyhedron
        delete: str
            delete polyhedra with this name
        list_polys: bool
            list existing polyhedra
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        return self._polyhedra(name, None, None, color, opacity, specular, specular_radius, ambient,
                               None, None, indices, delete, list_polys, all_frames)


    _parser_volume = ThrowingArgumentParser(prog="volume", description="read volumetric data from file")
    _parser_volume.add_argument("filename", nargs='?', help="File to read from. Text dap internal format, *.CHGCAR, *.PARCHG, *.WAVECAR, *.cube")
    _parser_volume.add_argument("-name", type=str, help="name to assign (default to filename)", default=None)
    _grp = _parser_volume.add_mutually_exclusive_group(required=True)
    _grp.add_argument("-delete", action='store', metavar="NAME", help="name of volume to delete", default=None)
    _grp.add_argument("-list_volumes", action='store_true', help="list existing volume representations")
    _grp.add_argument("-isosurface", type=float, action='store', metavar="THRESHOLD", help="isosurface threshold")
    _grp.add_argument("-volumetric", type=float, action='store', metavar="SCALE", help="volumetric value_to_opacity_factor")
    _parser_volume.add_argument("-color", nargs=3, type=float, metavar=("R", "G", "B"), help="color of representation", default=None)
    add_material_args_to_parser(_parser_volume, "representation")
    _parser_volume.add_argument("-file_args", nargs=argparse.REMAINDER, help="file-type specific arguments (use '-file_args -h' for file-type specific help)", default=[])
    #
    _volume_subparsers = {}
    #
    _volume_subparsers["CHGCAR"] = ThrowingArgumentParser(prog="file.CHGCAR -file_args", description="CHGCAR specific arguments")
    _volume_subparsers["CHGCAR"].add_argument("-component", choices=["total", "spin", "up", "down"], help="spin component to plot", default="total")
    #
    _volume_subparsers["WAVECAR"] = ThrowingArgumentParser(prog="WAVECAR", description="WAVECAR specific arguments")
    _volume_subparsers["WAVECAR"].add_argument("-nk", type=int, help="number of k-point (1-based)", required=True)
    _volume_subparsers["WAVECAR"].add_argument("-nband", type=int, help="number of band (1-based)", required=True)
    _volume_subparsers["WAVECAR"].add_argument("-component", choices=["up", "down"], help="spin component to plot, required for spin-polarized")
    _volume_subparsers["WAVECAR"].add_argument("-absolute_val", action='store_true', help="use absolute value rather than real-part of wavefunction")
    #
    _volume_subparsers["internal"] = ThrowingArgumentParser(prog="file.other -file_args", description="dap internal format specific arguments")
    _volume_subparsers["internal"].add_argument("-column", type=int, help="column (after indices) for dap internal text format", default=0)
    parsers["volume"] = _parser_volume
    def volume(self, filename=None, name=None, delete=None, list_volumes=None, isosurface=None, volumetric=None, color=None,
                   opacity=None, specular=None, specular_radius=None, ambient=None, file_args=None, _cmd_args=[]):
        """Visualize a volumetric dataset in periodic cell

        NEEDS_CMD_ARGS

        Parameters
        ----------
        filename: str
            file with database
        name: str
            name to associated with this volumetric rendering
        delete: str
            delete rendering with this name
        list_volumes: bool
            list existing renderings
        isosurface: float
            value at which to make isosurface
        volumetric: float
            coefficient for volute_to_opacity_factor
        color: (float, float, float)
            color of this rendering
        opacity: float
            isosurface rendering opacity
        specular: float
            isosurface rendering specular reflection
        specular_radius: float
            isosurface rendering specular shininess
        ambient: float
            isosurface rendering ambient light
        file_args: list(str):
            arguments to file-type-specific parser, in the form of strings
        """
        davtk_state = self.davtk_state

        if delete is not None:
            davtk_state.delete_volume_rep(davtk_state.cur_at(), delete)
        elif list_volumes:
            print("Defined volume names: "+" ".join(davtk_state.cur_at().volume_reps.keys()))
        else:
            creating_rep = isosurface is not None or volumetric is not None
            if filename is None and creating_rep:
                raise RuntimeError("volume got no filename")
            if name is None:
                if creating_rep:
                    name = filename
                else:
                    raise RuntimeError("not creating, deleting, or listing, need -name to modify")

            # create new property if needed
            if name not in davtk_state.volume_rep_prop:
                prop = new_prop( types.SimpleNamespace( color = (0.5, 0.5, 1.0), opacity = 0.5,
                    specular = 0.7, specular_radius = 0.1, ambient = 0.1 ) )
                davtk_state.volume_rep_prop[name] = prop

            update_prop(davtk_state.volume_rep_prop[name], types.SimpleNamespace(color=color, opacity=opacity,
                        specular=specular, specular_radius=specular_radius, ambient=ambient))

            if creating_rep:
                filename = Path(filename)
                if filename.name.endswith("CHGCAR") or filename.name.endswith("PARCHG"):
                    try:
                        sub_args = self._volume_subparsers["CHGCAR"].parse_args(file_args)
                    except ArgumentParserHelp:
                        return

                    chgcar = VaspChargeDensity(filename)
                    if sub_component == "total":
                        data = chgcar.chg[0]
                    else:
                        if not chgcar.is_spin_polarized():
                            raise RuntimeError("no spin-density available for non-spin_polarized calculation")
                        if sub_args.component == "spin":
                            data = chgcar.chgdiff[0]
                        elif sub_args.component == "up":
                            data = 0.5*(chgcar.chg[0] + chgcar.chgdiff[0])
                        elif sub_args.component == "down":
                            data = 0.5*(chgcar.chg[0] - chgcar.chgdiff[0])
                        else:
                            raise RuntimeError("volume CHGCAR should never get here")
                    data = np.ascontiguousarray(data.T)
                elif filename.name.endswith("WAVECAR"):
                    try:
                        sub_args = self._volume_subparsers["WAVECAR"].parse_args(file_args)
                    except ArgumentParserHelp:
                        return

                    from pymatgen.io.vasp.outputs import Wavecar
                    wf = Wavecar(filename) #, verbose=True)

                    if wf.spin == 2:
                        # spin polarized
                        if sub_args.component is None:
                            raise RuntimeError("component required for spin-polarized calculation")
                        if sub_args.component == "up":
                            wavecar = np.fft.ifftn(wf.fft_mesh(sub_args.nk - 1, sub_args.nband - 1, spin=0))
                        elif sub_args.component == "down":
                            wavecar = np.fft.ifftn(wf.fft_mesh(sub_args.nk - 1, sub_args.nband - 1, spin=1))
                        else:
                            assert "volume WAVECAR should never get here"
                    else:
                        # unpolarized
                        if sub_args.component is not None:
                            raise RuntimeError("component cannot be specified for unpolarized calculation")
                        wavecar = np.fft.ifftn(wf.fft_mesh(sub_args.nk - 1, sub_args.nband - 1))

                    if sub_args.absolute_val:
                        data = np.ascontiguousarray(np.abs(wavecar).T)
                    else:
                        data = np.ascontiguousarray(np.real(wavecar).T)

                    # normalize data so that corresponding density (vol. integral of square of abs) is normalized
                    norm_factor = np.product(data.shape) / davtk_state.cur_at().get_volume()
                    data *= np.sqrt(norm_factor) / np.linalg.norm(data)

                elif filename.name.endswith(".cube"):
                    data, _ = ase.io.cube.read_cube_data(filename)
                    data = np.ascontiguousarray(data.T)
                else:
                    try:
                        sub_args = self._volume_subparsers["internal"].parse_args(file_args)
                    except ArgumentParserHelp:
                        return

                    with open(filename) as fin:
                        extents = [int(i) for i in fin.readline().rstrip().split()]
                        if len(extents) != 3:
                            raise ValueError("Got bad number of extents {} != 3 on first line of '{}'".format(len(full_extents), filename))

                        # order of indices in data is being reversed
                        data = np.zeros(extents[::-1])
                        for l in fin:
                            fields = l.rstrip().split()
                            (i0, i1, i2) = fields[0:3]
                            v = fields[3 + sub_args.column]
                            try:
                                (i0, i1, i2) = (int(i0), int(i1), int(i2))
                            except:
                                (i0, i1, i2) = (int(np.round(float(i0)*extents[0])), int(np.round(float(i1)*extents[1])), int(np.round(float(i2)*extents[2])))
                            data[i2, i1, i0] = float(v)

                if isosurface is not None:
                    davtk_state.add_volume_rep(name, data, "isosurface", (isosurface,), ["volume"] + _cmd_args)
                if volumetric is not None:
                    raise RuntimeError("volume -volumetric not supported")
            # else: just modifying property

        return "cur"


    _parser_view = ThrowingArgumentParser(prog="view", description="set view position and orientation")
    _parser_view.add_argument("-lattice", help="Atoms.info key for lattice coordinates")
    _parser_view.add_argument("-direction", nargs='+', help="view and up directions with or without [] and () for Miller index notation")
    _parser_view.add_argument("-mag", type=float, help="view magnification (relative to current)", default=1.0)
    parsers["view"] = _parser_view
    def view(self, lattice=None, direction=None, mag=1.0):
        """Set view direction

        Parameters
        ----------
        lattice: str, optional
            Atoms.info key to use for Miller indices, otherwise Atoms.cell
        direction: list(str)
            two 3-vectors, first for view along direction, second for view up, in Cartesian
            (bare numbers) or Miller direction [] or plane normal ()
        mag: float
            magnification
        """
        if mag <= 0:
            raise ValueError("Can't set magnification <= 0")

        self.davtk_state.set_view(direction, lattice, mag)

        return "settings"


    _parser_alternate_cell_box = ThrowingArgumentParser(prog="alternate_cell_box",description="alternate (e.g. primitive) cell box to display")
    _parser_alternate_cell_box.add_argument("-all_frames", action="store_true",help="apply to all frames")
    _parser_alternate_cell_box.add_argument("-name","-n",help="name of info field", required=True)
    _grp = _parser_alternate_cell_box.add_mutually_exclusive_group(required=True)
    _grp.add_argument("-position","-p",type=float,nargs=3,help="Cartesian position of alternate cell box origin")
    _grp.add_argument("-atom","-a",type=int,help="Index of atom for alternate cell box origin")
    _grp.add_argument("-delete",action='store_true',help="Disable alternate cell box")
    parsers["alternate_cell_box"] = _parser_alternate_cell_box
    def alternate_cell_box(self, name, position=None, atom=None, delete=None, all_frames=False):
        """Display an alternate cell box (e.g. from a supercell or from a primitive cell

        Parameters
        ----------
        position: (float, float, float)
            cartesian position to offset box origin
        atom: int
            index of ato to offset box origin
        delete: bool
            delete alternate cell box
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        ats = self._get_ats(all_frames)

        for at in ats:
            if delete:
                if name in at.info["_vtk_alternate_cell_box"]:
                    del at.info["_vtk_alternate_cell_box"][name]
                else:
                    raise ValueError("-name {} not found".format(name))
            else:
                if atom is not None:
                    origin = atom
                else:
                    origin = position
            if name not in at.info:
                warnings.warn(f"alternate_cell_box info field {name} not in at.info")
            else:
                if "_vtk_alternate_cell_box" not in at.info:
                    at.info["_vtk_alternate_cell_box"] = {}
                at.info["_vtk_alternate_cell_box"][name] = origin

        return "cur"


    _parser_atom_override_type = ThrowingArgumentParser(prog="atom_override_type", description="override type of an atom")
    _parser_atom_override_type.add_argument("-all_frames", action="store_true", help="apply to all frames")
    _parser_atom_override_type.add_argument("-value", "-val", "-v", help="value to override with")
    _parser_atom_override_type.add_argument("-index", "-i", nargs='+', type=int, help="indices of atom(s) to override")
    _parser_atom_override_type.add_argument("-clear", "-c", action="store_true", help="clear overridden values for specified indices or all if not specified")
    parsers["atom_override_type"] = _parser_atom_override_type
    def atom_override_type(self, value=None, index=None, clear=False, all_frames=False):
        """Override the type of a specific atom

        Parameters
        ----------
        value: str
            value to override to
        index: list(int)
            indices of atoms whose type to override
        clear: bool
            clear overriding values
        all_frames: bool, default False
            apply to all frames rather than just current
        """
        if args.clear and args.value is not None:
            raise ValueError("Can't -clear and also set a value with -value")

        ats = self._get_ats(all_frames)

        for at in ats:
            if clear:
                if index is None:
                    del at.arrays["_vtk_override_type"]
                    continue
                else:
                    value="_NONE_"

            if "_vtk_override_type" not in at.arrays:
                at.new_array("_vtk_override_type", np.array(["_NONE_"] * len(at)).astype(object))
            at.arrays["_vtk_override_type"][np.array(index)] = value

        return "cur"
