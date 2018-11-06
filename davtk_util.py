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

class daVTK(object):
    def __init__(self, at_list, config):
        self.at_list = at_list
        self.config = config
        self.mappers = self.create_mappers()
        self.create_vtk_structures()
        self.cur_frame = 0

    def create_mappers(self):
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

    def update(self, frames=None):
        self.update_atoms(frames)

    def update_atoms(self, frames=None):

        for frame_i in self.frame_list(frames):
            at = self.at_list[frame_i]
            atom_type_array = get_atom_type_a(at)
            pos = at.get_positions()

            for i_at in range(len(at)):
                actor = at.arrays["_vtk_at_actor"][i_at]
                actor.SetMapper(self.mappers["sphere"])
                prop = get_atom_prop(self.config, atom_type_array[i_at], i_at, at.arrays)
                actor.SetProperty(prop)
                transform = vtk.vtkTransform()
                transform.Translate(pos[i_at])
                r = get_atom_radius(self.config, atom_type_array[i_at], i_at, at.arrays)
                transform.Scale(r, r, r)
                actor.SetUserMatrix(transform.GetMatrix())
                actor.i_at = i_at

    def create_vtk_structures(self):
        for at in self.at_list:
            at.arrays["_vtk_at_actor"] = [vtk.vtkActor() for i in range(len(at)) ]
            at.arrays["_vtk_picked"] = [False] * len(at)

        self.update_atoms()

    def set_shown_frame(self, renderer, dframe=None, frame_i=None):
        if dframe is not None:
            if frame_i is not None:
                raise ValueError("set_show_frame got both dframe and frame_i")
            self.cur_frame += dframe
        else:
            if frame_i is None:
                raise ValueError("set_show_frame got neither dframe and frame_i")

            self.cur_frame = frame_i

        # change list of shown actors
        renderer.RemoveAllViewProps()

        for actor in self.at_list[self.cur_frame].arrays["_vtk_at_actor"]:
            renderer.AddActor(actor)

        renderer.GetRenderWindow().Render()

