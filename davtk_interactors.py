import select, sys
import vtk

def toggle_pick(actor, picked_prop, at):
    if hasattr(actor,'prop_before_pick') and actor.prop_before_pick is not None:
        # restore to previous prop
        actor.SetProperty(actor.prop_before_pick)
        # reset saved prop to None
        actor.prop_before_pick = None
        at.arrays["_vtk_picked"][actor.i_at] = False
    else:
        # save current property, and set to picked prop
        actor.prop_before_pick = actor.GetProperty()
        actor.SetProperty(picked_prop)
        at.arrays["_vtk_picked"][actor.i_at] = True

class RubberbandSelect(vtk.vtkInteractorStyleRubberBand2D):
    def __init__(self,davtk_state,picked_prop,parent=None):
        self.AddObserver("LeftButtonReleaseEvent",self.leftButtonReleaseEvent)
        self.prev_style = None
        self.davtk_state = davtk_state
        self.picked_prop = picked_prop

    def set_prev_style(self, prev_style):
        self.prev_style = prev_style

    def leftButtonReleaseEvent(self,obj,event):
        p0 = self.GetStartPosition()
        p1 = self.GetEndPosition()

        picker = vtk.vtkAreaPicker()
        picker.AreaPick(p0[0], p0[1], p1[0], p1[1], self.GetDefaultRenderer())
        for p in picker.GetProp3Ds():
            toggle_pick(p, self.picked_prop, self.davtk_state.at_list[self.davtk_state.cur_frame])

        self.OnLeftButtonUp()

        if self.prev_style is None:
            raise ValueError("leftButtonReleaseEvent prev_style not set")

        self.GetInteractor().SetInteractorStyle(self.prev_style)
        self.prev_style.GetInteractor().Render()
        self.prev_Style = None

        return

class MouseInteractorHighLightActor(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self,config,davtk_state,select_style,parent=None):
        self.AddObserver("RightButtonPressEvent",self.rightButtonPressEvent)
        self.AddObserver("KeyPressEvent",self.keyPressEvent)
        self.AddObserver("TimerEvent",self.timerEvent)

        if(parent is not None):
            self.parent = parent
        else:
            self.parent = vtk.vtkRenderWindowInteractor()

        self.davtk_state = davtk_state
        self.config = config
        self.select_style = select_style

    def timerEvent(self,obj,event):
        if select.select([sys.stdin], [], [], 0)[0]:
            line = sys.stdin.readline()
            try:
                print "parsing line", line.rstrip()
                refresh = config_parse_line(self.config, line.rstrip())
                if refresh == "all":
                    self.davtk_state.update()
                elif refresh == "cur":
                    self.davtk_state.update("cur")
                elif refresh is not None:
                    raise ValueError("unknown refresh type "+str(refresh))
            except Exception, e:
                print "error parsing line",str(e)

            if self.GetInteractor() is not None:
                self.GetInteractor().Render()
            self.GetDefaultRenderer().GetRenderWindow().Render()

    def keyPressEvent(self,obj,event):
        k = self.parent.GetKeySym()
        if k == 's':
            self.GetInteractor().SetInteractorStyle(self.select_style)
            self.select_style.set_prev_style(self)
        elif k == 'd':
            self.davtk_state.delete_picked(self.GetDefaultRenderer())
        elif k == 'm':
            self.davtk_state.measure_picked()
        elif k == 'plus':
            self.davtk_state.set_shown_frame(renderer=self.GetDefaultRenderer(), dframe=1)
        elif k == 'minus':
            self.davtk_state.set_shown_frame(renderer=self.GetDefaultRenderer(), dframe=-1)
        elif k == 'c':
            if select.select([sys.stdin,],[],[],0.0)[0]:
                print "Have data!"
            else:
                print "No data"

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
            toggle_pick(self.NewPickedActor, self.select_style.picked_prop, self.davtk_state.at_list[self.davtk_state.cur_frame])
            if hasattr(self.NewPickedActor,'other_half'):
                toggle_pick(self.NewPickedActor.other_half, self.select_style.picked_prop, self.davtk_state.at_list[self.davtk_state.cur_frame])

        self.GetInteractor().Render()

        self.OnRightButtonDown()
        return
