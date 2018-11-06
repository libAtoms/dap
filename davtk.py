#!/usr/bin/env python

import ase.io, vtk, sys
from davtk_parse import config_parse_file
from davtk_util import *
from davtk_interactors import *
import argparse

cli_parse = argparse.ArgumentParser()
cli_parse.add_argument("--geometry",nargs=2,type=int,default=[800,800])
cli_parse.add_argument("files",nargs="+",type=str)
args = cli_parse.parse_args()

# read atoms from files
at_list = []
for f in args.files:
    for at in ase.io.read(f,":"):
        at.info["_vtk_filename"] = f
        at_list.append(at)

config = config_parse_file("davtk.config")

davtk_state = daVTK(at_list, config)

(min_pos, max_pos) = find_min_max(at_list)

# A renderer and render window
renderer = vtk.vtkRenderer()
renderer.SetBackground(config["background_color"])
renwin = vtk.vtkRenderWindow()
renwin.SetSize(args.geometry[0], args.geometry[1])
renwin.AddRenderer(renderer)

# An interactor for mouse stuff
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renwin)

# add the custom styles for regular interaction and area selection
sel_style = RubberbandSelect(davtk_state,config["picked_prop"],parent=interactor)
sel_style.SetDefaultRenderer(renderer)

def_style = MouseInteractorHighLightActor(config,davtk_state,sel_style,parent=interactor)
def_style.SetDefaultRenderer(renderer)
def_style.UseTimersOn()
interactor.CreateRepeatingTimer(100)

interactor.SetInteractorStyle(def_style)

# set up camera
camera = renderer.GetActiveCamera()
camera.ParallelProjectionOn()
camera.SetParallelScale(np.max(max_pos-min_pos))
transform = vtk.vtkTransform()
extent = max(max_pos) - min(min_pos)
camera.SetPosition([(max_pos[0]+min_pos[0])/2.0, (max_pos[1]+min_pos[1])/2.0, 1000.0+extent/2.0])
camera.SetFocalPoint((max_pos+min_pos)/2.0)
camera.SetClippingRange(1000-extent/2.0, 1000+3*extent/2.0)

# start viewing first frame
print "setting up first frame"
davtk_state.set_shown_frame(renderer, frame_i=0)

print "starting"
# Start
interactor.Initialize()
interactor.Start()
