from __future__ import print_function

import os, re, threading, queue
import vtk
from davtk.settings import DavTKSettings
from davtk.parse import parse_file
from davtk.state import *
from davtk.interactors import *

class Viewer(object):
    def __init__(self, at_list, win_size, win_name=None, init_commands=None):
        # read settings from home dir and current dir
        settings = DavTKSettings()
        try:
            parse_file(os.path.join(os.environ["HOME"],".daprc"), settings=settings)
        except IOError:
            pass
        try:
            parse_file(".daprc", settings=settings)
        except IOError:
            pass

        # A renderer and render window
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(settings["background_color"])
        renwin = vtk.vtkRenderWindow()
        renwin.SetSize(win_size[0], win_size[1])
        renwin.AddRenderer(renderer)

        # An interactor for mouse stuff
        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(renwin)

        self.davtk_state = DaVTKState(at_list, settings, renderer, self.interactor)

        # add the custom styles for regular interaction and area selection
        sel_style = RubberbandSelect(self.davtk_state,parent=self.interactor)
        sel_style.SetDefaultRenderer(renderer)

        def_style = MouseInteractorHighLightActor(settings,self.davtk_state,sel_style,parent=self.interactor)
        def_style.SetDefaultRenderer(renderer)
        def_style.UseTimersOn()
        self.interactor.CreateRepeatingTimer(100)

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
        camera.SetClippingRange(1000-extent/2.0, 1000+3*extent/2.0)

        l = vtk.vtkLight()
        l.SetPosition(100,70,200)
        l.SetLightTypeToCameraLight()
        renderer.GetLights().AddItem(l)

        # must do this rather late in the process
        renwin.SetWindowName(win_name)

        # read commands from atomic config headers
        for frame_i in range(len(self.davtk_state)):
            if "_vtk_commands" in at_list[frame_i].info:
                self.davtk_state.set_frame(str(frame_i))
                for cmd in at_list[frame_i].info["_vtk_commands"].split(";"):
                    parse_line(cmd, settings=settings, state=self.davtk_state)
                del at_list[frame_i].info["_vtk_commands"]

        self.davtk_state.set_frame("0")

        # now that atoms are read in and davtk_state exists, read any other commands (e.g. extra settings)
        for l in init_commands:
            for sub_l in l.split(";"):
                parse_line(sub_l, settings=settings, state=self.davtk_state)

        self.davtk_state.startup()

        self.davtk_state.prep_after_atoms_read()
        print("""DAP

        Use 'usage' for general usage info, and 'command -h' for detailed help on each command.
        Type 'h' in GUI window for GUI help message
        """)
        # Start
        self.interactor.Initialize()

        self.davtk_state.update("cur")

        self.cmd_queue = queue.Queue(0)
        self.davtk_state.cmd_queue = self.cmd_queue

    def get_cmd_queue(self):
        return self.cmd_queue

    def start(self, action, *args, **kwargs):

        self.cmd_thread = threading.Thread(target=action, args=args, kwargs=kwargs)
        self.cmd_thread.start()

        self.interactor.Start()
