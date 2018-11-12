import sys
import numpy as np
import vtk

def find_min_max(at_list):
    min_pos = sys.float_info.max
    max_pos = -sys.float_info.max

    for at in at_list:
        # update from positions
        p = at.get_positions()
        min_pos = np.minimum(min_pos, np.min(p, axis=0))
        max_pos = np.maximum(max_pos, np.max(p, axis=0))

        # update from cell extent
        c = at.get_cell()
        for i0 in range(2):
            for i1 in range(2):
                for i2 in range(2):
                    p = i0 * c[0] + i1 * c[1] + i2 * c[2]
                    min_pos = np.minimum(min_pos, p)
                    max_pos = np.maximum(max_pos, p)

    return (min_pos, max_pos)


def get_atom_type_a(at):
    if "atom_type" in at.arrays:
        atom_type = at.arrays["atom_type"]
    else:
        atom_type = [str(Z) for Z in at.get_atomic_numbers()]
    return atom_type

def get_atom_prop(settings, atom_type, i=None, arrays=None):
    if settings["atom_types"][atom_type]["colormap_func"] is not None and settings["atom_types"][atom_type]["colormap_field"] is not None:
        prop = vtk.vtkProperty()
        prop.DeepCopy(settings["atom_types"][atom_type]["prop"])
        prop.SetDiffuseColor(settings["atom_types"][atom_type]["colormap_func"](arrays[settings["atom_types"][atom_type]["colormap_field"]][i]))
        return prop
    else:
        return settings["atom_types"][atom_type]["prop"]

def get_atom_radius(settings, atom_type, i=None, arrays=None):
    if settings["atom_types"][atom_type]["radius_field"] is not None:
        return arrays[settings["atom_types"][atom_type]["radius_field"]][i]
    else:
        return settings["atom_types"][atom_type]["radius"]

class DaVTKState(object):
    def __init__(self, at_list, settings):
        self.at_list = at_list
        self.settings = settings
        self.mappers = self.create_shape_mappers()
        self.create_vtk_structures()
        self.cur_frame = 0

    def cur_at(self):
        return self.at_list[self.cur_frame]

    def create_shape_mappers(self):
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

    def frame_list(self, frames):
        if frames == "cur":
            return [self.cur_frame]
        elif isinstance(frames,int):
            return [frames]
        elif frames is None:
            return range(len(self.at_list))
        else:
            return frames

    def update(self, frames=None, atoms=None):
        self.update_atoms(frames, atoms)
        self.update_cell_boxes(frames)

    def update_atoms(self, frames=None, atoms=None):
        for frame_i in self.frame_list(frames):
            at = self.at_list[frame_i]
            atom_type_array = get_atom_type_a(at)
            pos = at.get_positions()

            if atoms is None:
                at_set = range(len(at))
            else:
                at_set = atoms

            for i_at in at_set:
                actor = at.arrays["_vtk_at_actor"][i_at]
                actor.SetMapper(self.mappers["sphere"])
                if at.arrays["_vtk_picked"][i_at]: 
                    prop = self.settings["picked_prop"]
                else:
                    prop = get_atom_prop(self.settings, atom_type_array[i_at], i_at, at.arrays)
                actor.SetProperty(prop)
                transform = vtk.vtkTransform()
                transform.Translate(pos[i_at])
                r = get_atom_radius(self.settings, atom_type_array[i_at], i_at, at.arrays)
                transform.Scale(r, r, r)
                actor.SetUserMatrix(transform.GetMatrix())

    def update_cell_boxes(self, frames=None):
        for frame_i in self.frame_list(frames):
            cell = self.at_list[frame_i].get_cell()

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

            actor = self.at_list[frame_i].info["_vtk_cell_box_actor"]
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(self.settings["cell_box_color"])
            actor.PickableOff()

    def create_vtk_structures(self):
        for at in self.at_list:
            at.arrays["_vtk_at_actor"] = [vtk.vtkActor() for i in range(len(at)) ]
            for (i_at, actor) in enumerate(at.arrays["_vtk_at_actor"]):
                actor.i_at = i_at
            at.arrays["_vtk_picked"] = np.array([False] * len(at))
            at.info["_vtk_cell_box_actor"] = vtk.vtkActor()

        self.update_atoms()
        self.update_cell_boxes()

    def set_shown_frame(self, renderer, dframe=None, frame_i=None):
        if dframe is not None:
            if frame_i is not None:
                raise ValueError("set_show_frame got both dframe and frame_i")
            self.cur_frame += dframe
        else: # dframe is None
            if frame_i is None:
                raise ValueError("set_show_frame got neither dframe and frame_i")
            self.cur_frame = frame_i

        # wrap around
        self.cur_frame = self.cur_frame % len(self.at_list)

        # change list of shown actors
        renderer.RemoveAllViewProps()

        # actors for atoms
        for actor in self.at_list[self.cur_frame].arrays["_vtk_at_actor"]:
            renderer.AddActor(actor)

        # actor for cell box
        renderer.AddActor(self.at_list[self.cur_frame].info["_vtk_cell_box_actor"])

        # need to do other actors, e.g. labels and bonds

        # refresh display
        renderer.GetRenderWindow().Render()

    def measure_picked(self):
        at = self.at_list[self.cur_frame]

        indices = np.where(at.arrays["_vtk_picked"])[0]
        p = at.get_positions()

        print "measure:"
        for i_at in indices:
            print "atom {} pos {} {} {}".format(i_at,p[i_at][0],p[i_at][1],p[i_at][2])
        for i in range(len(indices)):
            i_at = indices[i]
            for j_at in indices[i+1:]:
                dv = -at.get_distance(i_at,j_at,mic=True,vector=True)
                print "pair {} {} distance {} {} {} ({})".format(i_at,j_at,dv[0],dv[1],dv[2],np.linalg.norm(dv))

