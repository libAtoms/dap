#!/usr/bin/env python

import sys
import numpy as np
import vtk
import ase.io, ase.neighborlist
from davtk_parse import *

def create_mappers():
    mappers = {}

    source = vtk.vtkSphereSource()
    source.SetRadius(1.0)
    source.SetPhiResolution(8)
    source.SetThetaResolution(16)
    mappers["sphere"] = vtk.vtkPolyDataMapper()
    mappers["sphere"].SetInputConnection(source.GetOutputPort())

    source = vtk.vtkCylinderSource()
    source.SetRadius(1.0)
    source.SetHeight(1.0)
    source.SetResolution(8)
    mappers["cylinder"] = vtk.vtkPolyDataMapper()
    mappers["cylinder"].SetInputConnection(source.GetOutputPort())

    return mappers

def get_atom_type_a(at):
    if "atom_type" in at.arrays:
        atom_type = at.arrays["atom_type"]
    else:
        atom_type = [str(Z) for Z in at.get_atomic_numbers()]
    return atom_type

def get_atom_label(config, atom_type, i, arrays):
    if "atom_label" in arrays and len(arrays["atom_label"][i]) > 0:
        label = str(arrays["atom_label"][i])
    elif config["atom_types"][atom_type]["label_field"] is not None:
        if config["atom_types"][atom_type]["label_field"] == "NUMBER":
            label = str(i)
        else:
            label = str(arrays[config["atom_types"][atom_type]["label_field"]][i])
    else:
        label = str(i)
    return label

def get_atom_prop(config, atom_type, i=None, arrays=None):
    if config["atom_types"][atom_type]["colormap_func"] is not None and config["atom_types"][atom_type]["colormap_field"] is not None:
        prop = vtk.vtkProperty()
        prop.DeepCopy(config["atom_types"][atom_type]["prop"])
        prop.SetDiffuseColor(config["atom_types"][atom_type]["colormap_func"](arrays[config["atom_types"][atom_type]["colormap_field"]][i]))
        return prop
    else:
        return config["atom_types"][atom_type]["prop"]

def get_atom_radius(config, atom_type, i=None, arrays=None):
    if config["atom_types"][atom_type]["radius_field"] is not None:
        return arrays[config["atom_types"][atom_type]["radius_field"]][i]
    else:
        return config["atom_types"][atom_type]["radius"]

class Frame(object):
    def __init__(self):
        self.at = None
        self.at_actors = []
        self.bond_actors = []
        self.label_actors = []
        self.cell_box_actor = None
        return

    def measure_picked(self):
        print "measure:"
        indices = []
        for actor in self.at_actors:
            if hasattr(actor,'prop_before_pick') and actor.prop_before_pick is not None:
                indices.append(actor.i_at)
        p = self.at.get_positions()
        for i_at in indices:
            print "selected ",i_at, p[i_at]
        for i in range(len(indices)):
            i_at = indices[i]
            for j_at in indices[i+1:]:
                print "distance",i_at,j_at,self.at.get_distance(i_at,j_at,mic=True)

    def delete_atoms(self, renderer, indices):
        # delete actors
        for i_at in sorted(indices, reverse=True):
            renderer.RemoveActor(self.at_actors[i_at])
            del self.at_actors[i_at]
            renderer.RemoveActor(self.label_actors[i_at])
            del self.label_actors[i_at]
        # delete atoms
        del self.at[indices]
        # reset i_at components of atoms actors
        for i_at in range(len(self.at_actors)):
            self.at_actors[i_at].i_at = i_at

        renderer.GetRenderWindow().Render()

    def delete_picked(self, renderer):
        indices = []
        for actor in self.at_actors:
            if hasattr(actor,'prop_before_pick') and actor.prop_before_pick is not None:
                indices.append(actor.i_at)
        self.delete_atoms(renderer, indices)

    def activate(self, renderer):
        renderer.RemoveAllViewProps()

        if self.cell_box_actor is not None:
            renderer.AddActor(self.cell_box_actor)
        for (i, actor) in enumerate(self.at_actors):
            renderer.AddActor(actor)
        for actor in self.bond_actors:
            renderer.AddActor(actor)
        for actor in self.label_actors:
            renderer.AddActor(actor)

        renderer.GetRenderWindow().Render()

    def create_cell_box(self, cell, config):
        pts = vtk.vtkPoints()
        for i0 in range(2):
            for i1 in range(2):
                for i2 in range(2):
                    pts.InsertNextPoint(i0*cell[0] + i1*cell[1] + i2*cell[2])
        # 0 0 0    0 0 1    0 1 0    0 1 1     1 0 0  1 0 1   1 1 0   1 1 1
        lines = vtk.vtkCellArray()

        # origin to a1, a2, a3
        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,0)
        l.GetPointIds().SetId(1,1)
        lines.InsertNextCell(l)
        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,0)
        l.GetPointIds().SetId(1,2)
        lines.InsertNextCell(l)
        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,0)
        l.GetPointIds().SetId(1,4)
        lines.InsertNextCell(l)

        # a1, a2, a3 to 6 more
        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,1)
        l.GetPointIds().SetId(1,3)
        lines.InsertNextCell(l)
        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,1)
        l.GetPointIds().SetId(1,5)
        lines.InsertNextCell(l)

        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,2)
        l.GetPointIds().SetId(1,3)
        lines.InsertNextCell(l)
        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,2)
        l.GetPointIds().SetId(1,6)
        lines.InsertNextCell(l)

        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,4)
        l.GetPointIds().SetId(1,5)
        lines.InsertNextCell(l)
        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,4)
        l.GetPointIds().SetId(1,6)
        lines.InsertNextCell(l)

        # a1+a2+a3 to 3
        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,3)
        l.GetPointIds().SetId(1,7)
        lines.InsertNextCell(l)
        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,5)
        l.GetPointIds().SetId(1,7)
        lines.InsertNextCell(l)
        l = vtk.vtkLine()
        l.GetPointIds().SetId(0,6)
        l.GetPointIds().SetId(1,7)
        lines.InsertNextCell(l)

        linesPolyData = vtk.vtkPolyData()
        linesPolyData.SetPoints(pts)
        linesPolyData.SetLines(lines)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(linesPolyData)
        self.cell_box_actor = vtk.vtkActor()
        self.cell_box_actor.SetMapper(mapper)
        self.cell_box_actor.GetProperty().SetColor(config["cell_box_color"])
        self.cell_box_actor.PickableOff()

    def add_labels(self, at, config):
        self.label_actors = [None] * len(at)

        self.update_labels(config)

    def update_labels(self, config):
        atom_type_a = get_atom_type_a(self.at)
        pos = self.at.get_positions()
        for i_at in range(len(self.at)):
            label = vtk.vtkBillboardTextActor3D()
            label.SetInput(get_atom_label(config, atom_type_a[i_at], i_at, at.arrays))
            label.SetPosition(pos[i_at])
            r = get_atom_radius(config, atom_type_a[i_at], i_at, at.arrays)
            label.SetDisplayOffset(int(r*50),int(r*50))
            label.GetTextProperty().SetFontSize ( 24 )
            label.GetTextProperty().SetJustificationToLeft()
            self.label_actors[i_at] = label

    def add_atoms(self, at, config):
        # save ref to at
        self.at = at

        for i_at in range(len(at)):
            actor = vtk.vtkActor()
            self.at_actors.append(actor)

        self.update_atoms(config)

    def update_atoms(self, config):
        print "update_atoms"
        # get atom_type
        atom_type_a = get_atom_type_a(self.at)

        # update actors
        pos = self.at.get_positions()
        for i_at in range(len(at)):
            actor = self.at_actors[i_at]
            actor.SetMapper(mappers["sphere"])
            prop = get_atom_prop(config, atom_type_a[i_at], i_at, at.arrays)
            actor.SetProperty(prop)
            transform = vtk.vtkTransform()
            transform.Translate(pos[i_at])
            r = get_atom_radius(config, atom_type_a[i_at], i_at, at.arrays)
            transform.Scale(r, r, r)
            actor.SetUserMatrix(transform.GetMatrix())
            actor.i_at = i_at

def toggle_pick(actor):
    if hasattr(actor,'prop_before_pick') and actor.prop_before_pick is not None:
        # unpick
        actor.SetProperty(actor.prop_before_pick)
        actor.prop_before_pick = None
    else:
        # save previous property
        actor.prop_before_pick = actor.GetProperty()
        # create new picked property (should this be a previously allocated prop?)
        prop = vtk.vtkProperty()
        prop.SetColor(1.0, 1.0, 0.0)
        prop.SetDiffuse(1.0)
        prop.SetSpecular(0.0)
        actor.SetProperty(prop)

class RubberbandSelect(vtk.vtkInteractorStyleRubberBand2D):
    def __init__(self,parent=None):
        self.AddObserver("LeftButtonReleaseEvent",self.leftButtonReleaseEvent)
        self.prev_style = None

    def set_prev_style(self, prev_style):
        self.prev_style = prev_style

    def leftButtonReleaseEvent(self,obj,event):
        p0 = self.GetStartPosition()
        p1 = self.GetEndPosition()

        picker = vtk.vtkAreaPicker()
        picker.AreaPick(p0[0], p0[1], p1[0], p1[1], self.GetDefaultRenderer())
        for p in picker.GetProp3Ds():
            toggle_pick(p)

        self.OnLeftButtonUp()

        if self.prev_style is None:
            raise ValueError("leftButtonReleaseEvent prev_style not set")

        self.GetInteractor().SetInteractorStyle(self.prev_style)
        self.prev_style.GetInteractor().Render()
        self.prev_Style = None

        return

class MouseInteractorHighLightActor(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self,select_style,frames,parent=None):
        self.AddObserver("RightButtonPressEvent",self.rightButtonPressEvent)
        self.AddObserver("KeyPressEvent",self.keyPressEvent)

        if(parent is not None):
            self.parent = parent
        else:
            self.parent = vtk.vtkRenderWindowInteractor()

        self.select_style = select_style
        self.frames = frames
        self.cur_frame = 0


    def keyPressEvent(self,obj,event):
        k = self.parent.GetKeySym()
        if k == 's':
            self.GetInteractor().SetInteractorStyle(self.select_style)
            self.select_style.set_prev_style(self)
        elif k == 'd':
            self.frames[self.cur_frame].delete_picked(self.GetDefaultRenderer())
        elif k == 'm':
            self.frames[self.cur_frame].measure_picked()
        elif k == 'plus':
            self.cur_frame = (self.cur_frame+1) % len(self.frames)
            frames[self.cur_frame].activate(self.GetDefaultRenderer())
        elif k == 'minus':
            self.cur_frame = (self.cur_frame-1) % len(self.frames)
            frames[self.cur_frame].activate(self.GetDefaultRenderer())

        if self.GetInteractor() is not None:
            self.GetInteractor().Render()
        self.GetDefaultRenderer().GetRenderWindow().Render()

        # self.OnKeyPress()
        return

    def rightButtonPressEvent(self,obj,event):
        # get the actor at the picked position

        clickPos = self.GetInteractor().GetEventPosition()
        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())
        self.NewPickedActor = picker.GetActor()

        # If something was selected
        if self.NewPickedActor:
            toggle_pick(self.NewPickedActor)
            if hasattr(self.NewPickedActor,'other_half'):
                toggle_pick(self.NewPickedActor.other_half)

        self.GetInteractor().Render()

        self.OnRightButtonDown()
        return

# create mappers for spheres, cylinders, etc
mappers = create_mappers()

# read config
config = config_parse_file("davtk.config")
print config

# loop over configurations, adding them to frames, and keeping track of min, max pos for viewport
max_pos = np.array([-1.0e38] * 3)
min_pos = np.array([1.0e38] * 3)
frames = []
for f in sys.argv[1:]:
    this_file_ats = ase.io.read(f, ":")
    for at in this_file_ats:
        f = Frame()

        frames.append(f)

        f.add_atoms(at, config)
        f.add_labels(at, config)
        f.create_cell_box(at.get_cell(), config)

        # update min, max pos
        p = at.get_positions()
        c = at.get_cell()
        t_min_pos = np.min(p, axis=0)
        t_max_pos = np.max(p, axis=0)
        min_pos = np.minimum(min_pos, t_min_pos)
        max_pos = np.maximum(max_pos, t_max_pos)

        c = at.get_cell()
        for i0 in range(2):
            for i1 in range(2):
                for i2 in range(2):
                    p = i0 * c[0] + i1 * c[1] + i2 * c[2]
                    min_pos = np.minimum(min_pos, p)
                    max_pos = np.maximum(max_pos, p)


# A renderer and render window
renderer = vtk.vtkRenderer()
renderer.SetBackground(config["background_color"])
renwin = vtk.vtkRenderWindow()
renwin.SetSize(800,800)
renwin.AddRenderer(renderer)

# An interactor for mouse stuff
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renwin)

# add the custom styles for regular interaction and area selection
sel_style = RubberbandSelect(parent=interactor)
sel_style.SetDefaultRenderer(renderer)

def_style = MouseInteractorHighLightActor(parent=interactor,select_style=sel_style,frames=frames)
def_style.SetDefaultRenderer(renderer)

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
frames[0].activate(renderer)

print "starting"
# Start
interactor.Initialize()
interactor.Start()
