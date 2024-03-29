import sys, queue
import numpy as np
import vtk
from davtk.parse_utils import ArgumentParserHelp

from PyQt6.QtCore import QTimer

def pick_actors(at, actors, point_sets):
    new_bond_pick_statuses = {}
    for (actor, points) in zip(actors, point_sets):
        if hasattr(actor, "_vtk_type"):
            if actor._vtk_type in ["image_atom", "atoms_glyphs"]:
                if "_vtk_picked" not in at.arrays:
                    at.new_array("_vtk_picked",np.array([False]*len(at)))

                if actor._vtk_type == "image_atom":
                    i_at_list = [actor.i_at]
                elif actor._vtk_type == "atoms_glyphs":
                    i_at_list = [actor.i_at[point] for point in points]
                else:
                    raise RuntimeError("pick_actors() should never get here")

                at.arrays["_vtk_picked"][i_at_list] = np.logical_not(at.arrays["_vtk_picked"][i_at_list])
            elif actor._vtk_type == "bonds_glyphs":
                for point in points:
                    (i_at, i_bond) = actor.i_at_bond[point]
                    new_bond_pick_statuses[(i_at,i_bond)] = not at.bonds[i_at][i_bond]["picked"] 
            else:
                raise ValueError("picked something that's not an atom or a bond, rather: "+actor._vtk_type+"\n"+str(actor))
        else:
            raise ValueError("picked something that's not an atom or a bond, rather:\n"+str(actor))

    for ((i_at, i_bond), stat) in new_bond_pick_statuses.items():
        at.bonds.set_picked(i_at, i_bond, stat)

class RubberbandSelect(vtk.vtkInteractorStyleAreaSelectHover):
    def __init__(self, viewer, parent=None):
        self.AddObserver("LeftButtonReleaseEvent",self.leftButtonReleaseEvent)
        self.prev_style = None
        self.viewer = viewer

    def set_prev_style(self, prev_style):
        self.prev_style = prev_style

    def leftButtonReleaseEvent(self,obj,event):
        p0 = self.GetStartPosition()
        p1 = self.GetEndPosition()

        picker = vtk.vtkAreaPicker()
        picker.AreaPick(p0[0], p0[1], p1[0], p1[1], self.GetDefaultRenderer())
        at = self.viewer.davtk_state.cur_at()
        actors = picker.GetProp3Ds()
        points = []
        # find selected points
        for actor in actors:
            if not hasattr(actor, "point_to_input_point"):
                points.append([None])
                continue
            frustum = picker.GetFrustum()
            geom = vtk.vtkExtractGeometry()
            geom.SetImplicitFunction(frustum)
            geom.SetInputData(actor.GetMapper().GetInput())
            geom.Update()

            IDs = geom.GetOutputDataObject(0).GetPointData().GetArray("IDs")
            points.append([])
            for ID_i in range(IDs.GetNumberOfTuples()):
                points[-1].append(int(IDs.GetTuple(ID_i)[0]/actor.point_to_input_point))
            points[-1] = set(points[-1])

        pick_actors(self.viewer.davtk_state.cur_at(), actors, points)
        self.viewer.davtk_state.update()

        self.OnLeftButtonUp()

        if self.prev_style is None:
            raise ValueError("leftButtonReleaseEvent prev_style not set")

        self.GetInteractor().SetInteractorStyle(self.prev_style)
        self.prev_style.GetInteractor().Render()
        self.prev_Style = None

        return

class MouseInteractorHighLightActor(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self, viewer, select_style, parent=None):
        self.AddObserver("LeftButtonPressEvent",self.leftButtonPressEvent)
        self.AddObserver("LeftButtonReleaseEvent",self.leftButtonReleaseEvent)
        self.AddObserver("RightButtonPressEvent",self.rightButtonPressEvent)
        self.AddObserver("CharEvent",self.charEvent)
        # No need for observer when using QTimer instead (see note below)
        # self.AddObserver("TimerEvent",self.timerEvent)

        if parent is not None:
            self.parent = parent
        else:
            self.parent = vtk.vtkRenderWindowInteractor()

        # may slow things down?
        self.parent.AddObserver("ModifiedEvent",self.modifiedEvent)

        self.viewer = viewer
        self.select_style = select_style

        self.prev_size = (0,0)

        # timer via QTimer, as per
        #      https://stackoverflow.com/a/51608426
        # except that we _need_ to pass parent, apparently, otherwise
        # timer events aren't generated (or maybe just processed) except
        # on mouse click
        timer = QTimer(viewer.renwin)
        timer.timeout.connect(self.timerEvent)
        timer.start(100)


    def timerEvent(self, obj=None, event=None):
        try:
            line = self.viewer.davtk_state.cmd_queue.get(block=False)
        except queue.Empty:
            return
        try:
            refresh = self.viewer._parse_line(line)
            if refresh == "exit":
                sys.exit(0)
            else:
                self.viewer.davtk_state.update(what=refresh)
        except ArgumentParserHelp:
            pass
        except Exception as e:
            import traceback
            traceback.print_exc()
            print("error parsing line: ",e)

        if self.GetInteractor() is not None:
            self.GetInteractor().Render()
        self.GetDefaultRenderer().GetRenderWindow().Render()

    def modifiedEvent(self,obj,event):
        new_size = obj.GetSize()
        if new_size != self.prev_size:
            self.prev_size = new_size
            self.viewer.davtk_state.update()

    def charEvent(self,obj,event):
        k = self.parent.GetKeySym()
        if k == 'p':
            self.GetInteractor().SetInteractorStyle(self.select_style)
            self.select_style.set_prev_style(self)
        elif k == 'd':
            self.viewer.davtk_state.delete(self.viewer.cur_at(), atom_selection="picked", bond_selection="picked")
            self.viewer.davtk_state.update()
        elif k == 'm':
            self.viewer.davtk_state.measure()
        elif k == 'b':
            self.viewer.davtk_state.bond(self.viewer.cur_at(), name=None, at_type1='*', at_type2='*', criterion="picked")
            self.viewer.davtk_state.update()
        elif k == 'l':
            self.viewer.davtk_state.settings["atom_label"]["show"] = not self.viewer.davtk_state.settings["atom_label"]["show"]
            self.viewer.davtk_state.update()
        elif k == 'plus':
            self.viewer.davtk_state.update('+'+str(self.viewer.davtk_state.settings["frame_step"]))
        elif k == 'minus':
            self.viewer.davtk_state.update('-'+str(self.viewer.davtk_state.settings["frame_step"]))
        elif k in [ 'h', '?' ]:
            print("""GUI usage
h: help (this message)
+: next frame
-: prev frame
p: pick area mode (click and drag)
d: delete picked
m: measure picked
b: bond picked (only 2 atoms)
l: update labels
Mouse right click (two finger click on OS X): select point
Mouse left drag: trackball rotate
control-Mouse left drag: rotate in plane
shift-Mouse left drag: translate
Mouse scroll (two finger up/down drag on OS X): zoom
""")

        if self.GetInteractor() is not None:
            self.GetInteractor().Render()
        self.GetDefaultRenderer().GetRenderWindow().Render()

        # self.OnChar() # Forward Events to superclass handler

    def leftButtonPressEvent(self,obj,event):
        self.show_atom_labels_prev = self.viewer.davtk_state.settings["atom_label"]["show"]
        self.viewer.davtk_state.settings["atom_label"]["show"] = False

        self.show_legend_prev = self.viewer.davtk_state.settings["legend"]['show']
        self.viewer.davtk_state.settings["legend"]["show"] = False

        self.viewer.davtk_state.update("rotate")
        self.OnLeftButtonDown()
        return

    def leftButtonReleaseEvent(self,obj,event):
        self.viewer.davtk_state.settings["atom_label"]["show"] = self.show_atom_labels_prev

        self.viewer.davtk_state.settings["legend"]['show'] = self.show_legend_prev

        self.viewer.davtk_state.update("rotate")
        self.OnLeftButtonUp()

    def rightButtonPressEvent(self,obj,event):
        # get the actor at the picked position

        clickPos = self.GetInteractor().GetEventPosition()
        picker = vtk.vtkCellPicker()
        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())
        pickedActor = picker.GetActor()

        # If something was selected
        if pickedActor:
            at = self.viewer.davtk_state.cur_at()
            point = picker.GetPointId()
            try:
                point = int(point/pickedActor.point_to_input_point)
            except AttributeError:
                pass
            pick_actors(self.viewer.davtk_state.cur_at(), [pickedActor], [[point]])

        self.viewer.davtk_state.update()

        self.GetInteractor().Render()

        self.OnRightButtonDown()
        return
