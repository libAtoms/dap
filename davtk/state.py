import sys, ase.io, math
import numpy as np
import vtk
import ase.neighborlist
import re
from vtk.util.vtkImageImportFromArray import vtkImageImportFromArray
from vtk.util.numpy_support import vtk_to_numpy

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

def get_atom_type_a(settings, at):
    atom_type_field = settings["atom_type_field"]

    if atom_type_field == "Z":
        atom_type = [str(Z) for Z in at.get_atomic_numbers()]
    elif atom_type_field == "species":
        atom_type = [sp for sp in at.get_chemical_symbols()]
    else:
        atom_type = [str(val) for val in at.arrays[atom_type_field]]

    return atom_type

def get_atom_prop(settings, atom_type, i=None, arrays=None):
    if settings["atom_types"][atom_type]["colormap"] is not None:
        prop = vtk.vtkProperty()
        prop.DeepCopy(settings["atom_types"][atom_type]["prop"])
        colormap_field = settings["atom_types"][atom_type]["colormap"][1]
        colormap_func = settings["atom_types"][atom_type]["colormap"][2]
        prop.SetDiffuseColor(colormap_func(arrays[colormap_field][i]))
        return prop
    else:
        return settings["atom_types"][atom_type]["prop"]

def get_atom_radius(settings, atom_type, i=None, arrays=None):
    if settings["atom_types"][atom_type]["radius_field"] is not None:
        radius_field = settings["atom_types"][atom_type]["radius_field"][0]
        factor = settings["atom_types"][atom_type]["radius_field"][1]
        r = arrays[radius_field][i]*factor
    else:
        r = settings["atom_types"][atom_type]["radius"]
    if r is None:
        raise ValueError("Failed to find radius for atom type {}".format(atom_type))
    return r

class DavTKBonds(object):
    def __init__(self, at, settings):
        self.at = at
        self.settings = settings
        self.reinit()

    def __getitem__(self, key):
        return self.bonds[key]

    def reinit(self):
        self.bonds = [ [] for i in range(len(self.at)) ]

    def cutoff(self, name, in_cutoff, at_type, at_type2):

        def none_zero(x):
            return x if x is not None else 0.0

        def pair_match(at_type_a, i, j, at_type, at_type2):
            return ( ((at_type == '*' or str(at_type_a[i]) == at_type) and (at_type2 == '*' or str(at_type_a[j]) == at_type2)) or
                     ((at_type == '*' or str(at_type_a[j]) == at_type) and (at_type2 == '*' or str(at_type_a[i]) == at_type2)) )

        atom_type_array = get_atom_type_a(self.settings, self.at)

        if in_cutoff is None or len(in_cutoff) == 0: # fully auto
            max_cutoff = max([none_zero(self.settings["atom_types"][atom_type_array[i]]["bonding_radius"]) for i in range(len(self.at))])
            u_cutoff_min = lambda i1, i2 : 0.0
            u_cutoff_max = lambda i1, i2 : 0.5 * ( self.settings["atom_types"][atom_type_array[i1]]["bonding_radius"] +
                                                   self.settings["atom_types"][atom_type_array[i2]]["bonding_radius"] )
        elif len(in_cutoff) == 2: # min max
            max_cutoff = in_cutoff[1]
            u_cutoff_min = lambda i1, i2 : in_cutoff[0]
            u_cutoff_max = lambda i1, i2 : in_cutoff[1]
        elif len(in_cutoff) == 1: # max
            max_cutoff = in_cutoff[0]
            u_cutoff_min = lambda i1, i2 : 0.0
            u_cutoff_max = lambda i1, i2 : in_cutoff[0]

        if max_cutoff <= 0.0:
            return

        nn_list = ase.neighborlist.neighbor_list('ijDdS', self.at, max_cutoff, self_interaction=True)
        for (i, j, v, d, S) in zip(nn_list[0], nn_list[1], nn_list[2], nn_list[3], nn_list[4]):
            try:
                if d > 0.0 and d >= u_cutoff_min(i,j) and d <= u_cutoff_max(i, j) and pair_match(atom_type_array, i, j, at_type, at_type2):
                    self.bonds[i].append({ "j" : j, "v" : v, "d" : d, "S" : S, "name" : name, "picked" : False})
            except: # failed, presumably due to bonding_radius None
                pass

    def write_to_atoms_arrays(self, arrays_field="_vtk_bonds"):
        if arrays_field is None:
            arrays_field = "_vtk_bonds"
        bond_set_strs = []
        for (at_i, b_set) in enumerate(self.bonds):
            bond_l = []
            for b in b_set:
                v = [ b['j'], b['v'][0], b['v'][1], b['v'][2], b['S'][0], b['S'][1], b['S'][2], b['picked'], b['name'] ]
                bond_l.append("_".join([str(vi) for vi in v]))
            bond_set_strs.append(",".join(bond_l) if len(bond_l) > 0 else "_NONE_")
        self.at.arrays[arrays_field] = np.array(bond_set_strs)

    def read_from_atoms_arrays(self, arrays_field="_vtk_bonds"):
        if arrays_field is None:
            arrays_field = "_vtk_bonds"
        self.reinit()
        def str_to_bool(s):
            if s in ["F", "False", "false", "f"]:
                return False
            elif s in ["T", "True", "true", "t"]:
                return True
            else:
                raise ValueError("str_to_bool cannot parse '{}'".format(s))
        for (at_i, b_str) in enumerate(self.at.arrays[arrays_field]):
            if b_str == "_NONE_":
                continue
            bond_strs = b_str.split(",")
            for bond_str in bond_strs:
                m = re.search('^([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_(.*)$', bond_str)
                if m:
                    (j, v, S, picked, name) = ( int(m.group(1)), [float(m.group(i)) for i in range(2,5)], 
                                                [int(m.group(i)) for i in range(5,8)], str_to_bool(m.group(8)), m.group(9) )
                    if name not in self.settings["bond_types"]:
                        raise ValueError("Unknown bond_type '{}' reading from Atoms object, known types: {}\n".format(name, list(self.settings["bond_types"].keys()))+
                                         "            Perhaps you forgot to restore saved settings?")
                    self.bonds[at_i].append( {"j" : j, "v" : v, "d" : np.linalg.norm(v), "S" : S, "name" : name, "picked" : picked } )

    def pair_mic(self, name, ind1, ind2):
        v = self.at.get_distance(ind1, ind2, mic=True, vector=True)
        self.bonds[ind1].append( {"j" : ind2, "v" : v, "d" : np.linalg.norm(v), "S" : [0,0,0], "name" : name, "picked" : False } )
        self.bonds[ind2].append( {"j" : ind1, "v" : -v, "d" : np.linalg.norm(v), "S" : [0,0,0], "name" : name, "picked" : False } )

    def set_picked(self, i_at, j_ind, stat):
        b = self.bonds[i_at][j_ind]
        b["picked"] = stat
        j_at = b["j"]
        if j_at != i_at:
            for (bb_i, bb) in enumerate(self.bonds[j_at]):
                if bb["j"] == i_at and np.all(b["S"] == -bb["S"]):
                    self.bonds[b["j"]][bb_i]["picked"] = stat
                    return
        raise ValueError("set_picked failed to find opposite for {} {}".format(i_at, j_at))

    def delete_one(self, i_at, j_ind):
        b = self.bonds[i_at][j_ind]
        del self.bonds[i_at][j_ind]
        j_at = b["j"]
        if j_at != i_at:
            for (bb_i, bb) in enumerate(self.bonds[j_at]):
                if bb["j"] == i_at and np.all(b["S"] == -bb["S"]):
                    del self.bonds[b["j"]][bb_i]
                    return
        raise ValueError("delete_one failed to find opposite for {} {}".format(i_at, j_at))

    def delete_atoms(self, at_inds):
        orig_n = len(self.bonds)

        new_indices = np.array( [-1] * orig_n )
        new_indices[self.at.arrays["_vtk_orig_indices"]] = range(len(self.at))

        # remove records for atoms
        for i_at in sorted(at_inds, reverse=True):
            del self.bonds[i_at]
        # remove records that refer to atoms
        for i_at in range(len(self.at)):
            del_i_bonds = []
            for (i_bond, b) in enumerate(self.bonds[i_at]):
                if new_indices[b["j"]] < 0:
                    del_i_bonds.append(i_bond)
                else:
                    b["j"] = new_indices[b["j"]]
            for j in sorted(del_i_bonds, reverse=True):
                del self.bonds[i_at][j]

class DaVTKState(object):
    def __init__(self, at_list, settings, renderer, iRen=None):
        self.at_list = at_list
        self.settings = settings
        self.mappers = self.create_shape_mappers()
        self.create_vtk_structures()
        self.cur_frame = 0
        self.renderer = renderer
        self.iRen = iRen

        self.atom_actor_pool = []
        self.cur_n_atom_actors = 0
        self.label_actor_pool = []
        self.cur_n_label_actors = 0
        self.bond_actor_pool = []
        self.cur_n_bond_actors = 0

        self.saved_views = {}

        self.update()

    def __len__(self):
        return len(self.at_list)

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
        elif frames is None or frames == "all":
            return range(len(self.at_list))
        else:
            return frames

    def delete(self, atoms=None, bonds=None, frames="cur"):
        for frame_i in self.frame_list(frames):

            # delete atoms and their bonds
            at = self.at_list[frame_i]
            at.arrays["_vtk_orig_indices"] = np.array(range(len(at)))
            if atoms is not None:
                if atoms == "picked":
                    at_inds = np.where(at.arrays["_vtk_picked"])[0]
                elif isinstance(atoms,int):
                    at_inds = np.array([atoms])
                else:
                    at_inds = np.array(atoms)
                del at[at_inds]
                if hasattr(at, "bonds"):
                    at.bonds.delete_atoms(at_inds)

            # delete requested bonds
            if bonds is not None and hasattr(at, "bonds"):
                if bonds == "picked":
                    for i_at in range(len(at)):
                        for j_ind in [jj for jj in range(len(at.bonds[i_at])) if at.bonds[i_at][jj]["picked"]]:
                            at.bonds.delete_one(i_at, j_ind)
                elif bonds == "all":
                    at.bonds.reinit()
                else:
                    raise ValueError("delete bonds only accepts 'picked' or 'all' or 'all'")

        self.update(frames)
        self.show_frame(dframe=0)

    def update(self, frames=None):
        self.update_cell_boxes(frames)
        self.show_frame(dframe=0)

    def update_bonds(self, at):
        if not hasattr(at, "bonds"):
            for i in range(self.cur_n_bond_actors):
                self.bond_actor_pool[i].SetVisibility(False)
            self.cur_n_bond_actors = 0
            return

        n_bonds = 0
        for i_at in range(len(at)):
            for (i_bond, b) in enumerate(at.bonds[i_at]):
                if i_at < b["j"]:
                    continue
                n_bonds += 1

        if 2*n_bonds > len(self.bond_actor_pool):
            prev_pool_size = len(self.bond_actor_pool)
            self.bond_actor_pool.extend([vtk.vtkActor() for i in range(2*n_bonds-prev_pool_size)])
            for actor in self.bond_actor_pool[prev_pool_size:]:
                actor.SetMapper(self.mappers["cylinder"])

        pos = at.get_positions()

        global_i_bond = 0
        for i_at in range(len(at)):
            for (i_bond, b) in enumerate(at.bonds[i_at]):
                i = i_at
                j = b["j"]
                if i < j:
                    continue
                dr = np.array(b["v"])
                dr_norm = b["d"]
                picked = b["picked"]
                name = b["name"]

                rad = self.settings["bond_types"][name]["radius"]

                axis = np.cross(dr, [0.0, 1.0, 0.0])
                if np.linalg.norm(axis) == 0.0:
                    axis = None
                else:
                    axis /= np.linalg.norm(axis)
                dr_hat = dr / dr_norm
                angle = -np.arccos(np.dot(dr_hat, [0.0, 1.0, 0.0]))*180.0/np.pi

                # first half
                actor_1 = self.bond_actor_pool[global_i_bond]

                transform = vtk.vtkTransform()
                transform.Translate(pos[i]+dr/4.0)
                if axis is not None:
                    transform.RotateWXYZ(angle, axis)
                transform.Scale(rad, dr_norm/2.0, rad)
                actor_1.SetUserMatrix(transform.GetMatrix())

                # second half bond
                actor_2 = self.bond_actor_pool[global_i_bond+1]

                transform = vtk.vtkTransform()
                transform.Translate(pos[j]-dr/4.0)
                if axis is not None:
                    transform.RotateWXYZ(angle, axis)
                transform.Scale(rad, dr_norm/2.0, rad)
                actor_2.SetUserMatrix(transform.GetMatrix())


                if picked:
                    actor_1.SetProperty(self.settings["picked_prop"])
                    actor_2.SetProperty(self.settings["picked_prop"])
                else:
                    actor_1.SetProperty(self.settings["bond_types"][name]["prop"])
                    actor_2.SetProperty(self.settings["bond_types"][name]["prop"])

                actor_1.i_at_bond = (i_at, i_bond)
                actor_2.i_at_bond = (i_at, i_bond)

                actor_1.SetVisibility(True)
                actor_2.SetVisibility(True)

                global_i_bond += 2

        for i in range(n_bonds*2, self.cur_n_bond_actors):
            self.bond_actor_pool[i].SetVisibility(False)

        self.cur_n_bond_actors = 2*n_bonds

    def label_offset_world(self, pos):
        self.renderer.GetActiveCamera()

        # screen pos of arbitrary world point
        pt_disp = [0.0, 0.0, 0.0]
        self.iRen.GetInteractorStyle().ComputeWorldToDisplay(self.renderer, pos[0], pos[1], pos[2], pt_disp)
        # world pos of point 10 pixel away
        pt_disp[0] += 10.0
        pt_world = [0.0, 0.0, 0.0, 0.0]
        self.iRen.GetInteractorStyle().ComputeDisplayToWorld(self.renderer, pt_disp[0], pt_disp[1], pt_disp[2], pt_world)
        dp_world = np.linalg.norm(np.array(pos)-np.array(pt_world[0:3]))/10

        return dp_world

    def update_labels(self, at):
        if len(at) > len(self.label_actor_pool):
            prev_pool_size = len(self.label_actor_pool)
            self.label_actor_pool.extend([vtk.vtkBillboardTextActor3D() for i in range(len(at)-prev_pool_size)])

        atom_type_array = get_atom_type_a(self.settings, at)
        pos = at.get_positions()

        dp_world = self.label_offset_world(pos[0])

        for i_at in range(len(at)):
            label_actor = self.label_actor_pool[i_at]
            label_str = None
            if "_vtk_label" in at.arrays:
                if at.arrays["_vtk_label"][i_at] == "_NONE_":
                    label_str = ""
                elif len(at.arrays["_vtk_label"][i_at]) == 0 or at.arrays["_vtk_label"][i_at] == "'''" or at.arrays["_vtk_label"][i_at] == '""':
                    label_str = None
                else:
                    label_str = at.arrays["_vtk_label"][i_at]
            if label_str is None:
                label_field = self.settings["atom_types"][atom_type_array[i_at]]["label"]
                if label_field is None or label_field == "ID":
                    label_str = str(i_at)
                elif label_field == "Z" or label_field == "atomic_number":
                    label_str = str(at.get_atomic_numbers()[i_at])
                elif label_field == "species":
                    label_str = at.get_chemical_symbols()[i_at]
                else:
                    label_str = str(at.arrays[label_field][i_at])
            label_actor.SetInput(label_str)
            label_actor.SetPosition(pos[i_at])
            r = get_atom_radius(self.settings, atom_type_array[i_at], i_at, at.arrays)
            if dp_world > 0:
                dp_disp = 0.7*r/dp_world
            else:
                dp_disp = 0
            label_actor.SetDisplayOffset(int(dp_disp), int(dp_disp))
            label_actor.SetTextProperty(self.settings["label_text_prop"])
            label_actor.PickableOff()
            label_actor.SetVisibility(True)

        for i in range(len(at), self.cur_n_label_actors):
            self.label_actor_pool[i].SetVisibility(False)
        self.cur_n_label_actors = len(at)

    def update_atoms(self, at):
        if len(at) > len(self.atom_actor_pool):
            prev_pool_size = len(self.atom_actor_pool)
            self.atom_actor_pool.extend([vtk.vtkActor() for i in range(len(at)-prev_pool_size)])
            for actor in self.atom_actor_pool[prev_pool_size:]:
                actor.SetMapper(self.mappers["sphere"])

        atom_type_array = get_atom_type_a(self.settings, at)
        pos = at.get_positions()

        for i_at in range(len(at)):
            actor = self.atom_actor_pool[i_at]
            if at.arrays["_vtk_picked"][i_at]: 
                prop = self.settings["picked_prop"]
            else:
                prop = get_atom_prop(self.settings, atom_type_array[i_at], i_at, at.arrays)
            actor.SetProperty(prop)
            transform = vtk.vtkTransform()
            transform.Translate(pos[i_at])
            r = get_atom_radius(self.settings, atom_type_array[i_at], i_at, at.arrays)
            actor.r = r
            transform.Scale(r,r,r)
            actor.SetUserMatrix(transform.GetMatrix())
            # update in case numbers changed
            actor.i_at = i_at
            actor.SetVisibility(True)

        for i in range(len(at), self.cur_n_atom_actors):
            self.atom_actor_pool[i].SetVisibility(False)
        self.cur_n_atom_actors = len(at)

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

            actor = self.at_list[frame_i]._NOPRINT_vtk_cell_box_actor
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(self.settings["cell_box_color"])
            actor.PickableOff()

    def create_vtk_structures(self, frames=None):
        for frame_i in self.frame_list(frames):
            at = self.at_list[frame_i]

            at.arrays["_vtk_picked"] = np.array([False] * len(at))
            at._NOPRINT_vtk_cell_box_actor = vtk.vtkActor()
            at.info["_vtk_show_labels"] = False

    def show_image_atoms(self, at, pos):
        cell = at.get_cell()
        cell_inv = at.get_reciprocal_cell().T
        for (i_at, actor) in enumerate(self.atom_actor_pool[0:self.cur_n_atom_actors]):
            # real image
            self.renderer.AddActor(actor)
            # periodic images if needed
            if "_vtk_images" in at.info:
                n0 = int(math.ceil(at.info["_vtk_images"][0]))
                n1 = int(math.ceil(at.info["_vtk_images"][1]))
                n2 = int(math.ceil(at.info["_vtk_images"][2]))
                for i0 in range(-n0, n0+1):
                    for i1 in range(-n1, n1+1):
                        for i2 in range(-n2, n2+1):
                            if (i0,i1,i2) == (0,0,0):
                                continue
                            img_pos = pos[i_at] + np.dot([i0, i1, i2], cell)
                            img_pos_scaled = np.dot(img_pos, cell_inv)
                            if (img_pos_scaled[0] >= -at.info["_vtk_images"][0] and
                                img_pos_scaled[0] <= 1+at.info["_vtk_images"][0] and
                                img_pos_scaled[1] >= -at.info["_vtk_images"][1] and
                                img_pos_scaled[1] <= 1+at.info["_vtk_images"][1] and
                                img_pos_scaled[2] >= -at.info["_vtk_images"][2] and
                                img_pos_scaled[2] <= 1+at.info["_vtk_images"][2]):
                               img_actor = vtk.vtkActor()
                               img_actor.SetProperty(actor.GetProperty())
                               img_actor.SetMapper(actor.GetMapper())
                               transform = vtk.vtkTransform()
                               transform.Translate(img_pos)
                               transform.Scale(actor.r, actor.r, actor.r)
                               img_actor.SetUserMatrix(transform.GetMatrix())
                               img_actor.i_at = i_at
                               self.renderer.AddActor(img_actor)


    def show_legend(self, at, pos):
        display_size = self.renderer.GetRenderWindow().GetSize()
        dp_world = self.label_offset_world(pos[0])
        unique_atom_types = sorted(list(set(get_atom_type_a(self.settings, at))))
        legend_sphere_actors = [ vtk.vtkActor() for i in range(len(unique_atom_types)) ]
        legend_label_actors = [ vtk.vtkActor() for i in range(len(unique_atom_types)) ]

        for (l_i, (legend_sphere_actor, legend_label_actor, atom_type)) in enumerate(zip(legend_sphere_actors, legend_label_actors, unique_atom_types)):
            atom_i_of_type = unique_atom_types.index(atom_type)

            legend_sphere_actor.SetMapper(self.mappers["sphere"])
            prop = get_atom_prop(self.settings, atom_type, atom_i_of_type, at.arrays)
            legend_sphere_actor.SetProperty(prop)
            legend_sphere_actor.SetScale(self.atom_actor_pool[atom_i_of_type].r,
                                         self.atom_actor_pool[atom_i_of_type].r,
                                         self.atom_actor_pool[atom_i_of_type].r)
            legend_pos = self.settings["legend"]['position'] % display_size

            coord = vtk.vtkCoordinate()
            coord.SetCoordinateSystemToDisplay()
            coord.SetValue(legend_pos[0], legend_pos[1]-self.settings["legend"]['spacing']*l_i, 0.01)
            sphere_pos_world = coord.GetComputedWorldValue(self.renderer)

            legend_sphere_actor.SetPosition(sphere_pos_world[0:3])
            legend_sphere_actor.PickableOff()
            legend_sphere_actor.VisibilityOn()
            self.renderer.AddActor(legend_sphere_actor)

            legend_label_actor = vtk.vtkBillboardTextActor3D()
            legend_label_actor.SetInput(atom_type)
            legend_label_actor.SetPosition(sphere_pos_world[0:3])
            if dp_world > 0:
                dp_disp = self.atom_actor_pool[atom_i_of_type].r/dp_world
            else:
                dp_disp = 0
            legend_label_actor.SetDisplayOffset(int(1.3*dp_disp), -int(dp_disp/2.0))
            legend_label_actor.SetTextProperty(self.settings["label_text_prop"])
            legend_label_actor.PickableOff()
            legend_label_actor.VisibilityOn()
            self.renderer.AddActor(legend_label_actor)

    def set_frame(self, dframe=None, frame_i=None):
        if dframe is not None:
            if frame_i is not None:
                raise ValueError("set_show_frame got both dframe and frame_i")
            self.cur_frame += dframe
        else: # dframe is None
            if frame_i is None:
                raise ValueError("set_show_frame got neither dframe and frame_i")
            if frame_i < 0 or frame_i >= len(self.at_list):
                raise ValueError("set_show_frame got frame_i {} out of range {} --- {}".format(frame_i, 0, len(self.at_list)))
            self.cur_frame = frame_i

        self.renderer.SetBackground(self.settings["background_color"])

        # wrap around
        self.cur_frame = self.cur_frame % len(self.at_list)

    def show_frame(self, dframe=None, frame_i=None):
        self.set_frame(dframe, frame_i)

        at = self.at_list[self.cur_frame]

        self.update_atoms(at)
        self.update_labels(at)
        self.update_bonds(at)

        # remove all existing actors
        self.renderer.RemoveAllViewProps()

        # actor for cell box
        self.renderer.AddActor(at._NOPRINT_vtk_cell_box_actor)

        # create actors for atom images
        pos = at.get_positions()
        self.show_image_atoms(at, pos)

        # need to do other actors, e.g. labels and bonds
        for actor in self.bond_actor_pool[0:self.cur_n_bond_actors]:
            self.renderer.AddActor(actor)

        if self.settings["legend"]['show']:
            self.show_legend(at, pos)

        if at.info.get("_vtk_show_labels", False):
            for actor in self.label_actor_pool[0:self.cur_n_label_actors]:
                self.renderer.AddActor(actor)

        try:
            for volume_rep in at.volume_reps.values():
                for (actor, _) in volume_rep[2]:
                    self.renderer.AddActor(actor)
        except AttributeError:
            pass

        # config_n
        config_n_actor = vtk.vtkTextActor()
        config_n_actor.SetInput(str(self.cur_frame))
        config_n_actor.SetTextProperty(self.settings["config_n_text_prop"])
        config_n_actor.SetDisplayPosition(20,20)
        self.renderer.AddActor(config_n_actor)

        # refresh display
        self.renderer.GetRenderWindow().Render()

    def measure(self, n=None, frame_i=None):
        if frame_i is None:
            at = self.at_list[self.cur_frame]
        else:
            at = self.at_list[frame_i]

        if n is None:
            at_indices = np.where(at.arrays["_vtk_picked"])[0]
        else:
            at_indices = n

        p = at.get_positions()

        print("measure:")
        for i_at in at_indices:
            print("atom {} pos {} {} {}".format(i_at,p[i_at][0],p[i_at][1],p[i_at][2]))
        for i in range(len(at_indices)):
            i_at = at_indices[i]
            for j_at in at_indices[i+1:]:
                dv = -at.get_distance(i_at,j_at,mic=True,vector=True)
                print("atom-distance {} {} vec {} {} {} ({})".format(i_at,j_at,dv[0],dv[1],dv[2],np.linalg.norm(dv)))

        if n is None and hasattr(at, "bonds"):
            for i_at in range(len(at)):
                for b in at.bonds[i_at]:
                    if b["picked"] and b["j"] >= i_at:
                        print("bond {} {} vec {} {} {} ({})".format(i_at, b["j"], b["v"][0], b["v"][1], b["v"][2], b["d"]))
            for i_at in range(len(at)):
                for b in at.bonds[i_at]:
                    j_at = b["j"]
                    if b["picked"]:
                        for bb in at.bonds[j_at]:
                            k_at = bb["j"]
                            if bb["picked"] and k_at <= i_at and (k_at != i_at or np.any(b["S"] != -bb["S"])):
                                ang = 180.0/np.pi * np.arccos(np.dot(b["v"],-bb["v"])/(np.linalg.norm(b["v"])*np.linalg.norm(bb["v"])))
                                print("bond-angle {} {} {} angle {}".format(i_at, j_at, k_at, ang))

    def duplicate(self, n_dup, frames=None):
        for frame_i in self.frame_list(frames):
            at = self.at_list[frame_i]
            # store original info if this is first duplicate
            if "dup_orig_n" not in at.info:
                at.info["dup_orig_n"] = len(at)
            if "dup_orig_cell" not in at.info:
                at.info["dup_orig_cell"] = at.get_cell()
            if "dup_orig_index" not in at.arrays:
                at.new_array("dup_orig_index",np.array(range(len(at))))
            # store current info and duplicate atoms
            prev_n = len(at)
            prev_cell = at.get_cell()
            prev_pos = at.get_positions()
            p = list(prev_pos)
            for i0 in range(n_dup[0]):
                for i1 in range(n_dup[1]):
                    for i2 in range(n_dup[2]):
                        if (i0, i1, i2) == (0, 0, 0):
                            continue
                        at.extend(at[0:prev_n])
                        p.extend(prev_pos + np.dot([i0,i1,i2], prev_cell))
            at.set_positions(p)
            at.set_cell(np.dot(np.diagflat([n_dup[0], n_dup[1], n_dup[2]]), prev_cell), False)
            # need to duplicate bonds
            if hasattr(at, "bonds"):
                at.bonds.reinit()

        self.create_vtk_structures(frames)
        self.update(frames)

    def bond(self, name, at_type1, at_type2, criterion, frames=None):
        if name is None:
            name = self.settings["default_bond_type"]

        for frame_i in self.frame_list(frames):
            at = self.at_list[frame_i]

            if not hasattr(at, "bonds"):
                at.bonds = DavTKBonds(at, self.settings)

            if criterion == "picked":
                indices = np.where(self.cur_at().arrays["_vtk_picked"])[0]
                if len(indices) != 2:
                    raise ValueError("tried to bond picked, but number of atoms picked {} != 2".format(len(indices)))
                at.bonds.pair_mic(name, indices[0], indices[1])
                at.arrays["_vtk_picked"][indices] = False
            elif criterion[0] == "n":
                indices = criterion[1]
                if len(indices) != 2:
                    raise ValueError("tried to bond pair, but number of indices passed {} != 2".format(len(indices)))
                at.bonds.pair_mic(name, indices[0], indices[1])
            elif criterion[0] == "cutoff":
                at.bonds.cutoff(name, criterion[1], at_type1, at_type2)
            else:
                raise ValueError("Unknown bonding criterion type '{}'".format(criterion[0]))

        self.update(frames)

    def snapshot(self, filename=None, mag=1):
        renderLarge = vtk.vtkRenderLargeImage()
        renderLarge.SetInput(self.renderer)
        renderLarge.SetMagnification(mag)

        if filename is not None:
            writer = vtk.vtkPNGWriter()
            writer.SetInputConnection(renderLarge.GetOutputPort())
            writer.SetFileName(filename)
            writer.Write()
            return None
        else:
            renderLarge.Update()
            img = renderLarge.GetOutput()
            (width, height, _) = img.GetDimensions()
            raw_data = img.GetPointData().GetScalars()
            return vtk_to_numpy(raw_data).reshape(height,width,3)

    def get_view(self):
        cam = self.renderer.GetActiveCamera()
        quantities = [cam.GetParallelScale()] + list(cam.GetPosition()) + list(cam.GetFocalPoint()) + list(cam.GetClippingRange()) + list(cam.GetViewUp())
        return quantities

    def restore_view(self, quantities):
        if len(quantities) != 12:
            raise ValueError("restore_view got wrong number of values {} != 12".format(len(quantities)))
        cam = self.renderer.GetActiveCamera()
        cam.SetParallelScale(quantities[0])
        cam.SetPosition(quantities[1:4])
        cam.SetFocalPoint(quantities[4:7])
        cam.SetClippingRange(quantities[7:9])
        cam.SetViewUp(quantities[9:12])

    def startup(self):
        # restore view

        self.show_frame(frame_i=0)

    def array_to_image(self, data):
        img = vtkImageImportFromArray()
        img.SetArray(data)
        img.Update()
        return img

    def delete_volume_rep(self, name):
        del self.cur_at().volume_reps[name]

    def make_isosurface(self, img, transform, params):

        isosurface = vtk.vtkMarchingCubes()
        isosurface.SetInputData( img.GetOutput() )
        isosurface.ComputeNormalsOn()
        isosurface.SetValue( 0, params[0] )
        isosurface.Update()

        scaled_iso = vtk.vtkTransformFilter()
        scaled_iso.SetInputData( isosurface.GetOutput() )

        scaled_iso.SetTransform(transform)
        scaled_iso.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData( scaled_iso.GetOutput() )
        mapper.ScalarVisibilityOff()
        mapper.Update()

        actor = vtk.vtkActor()
        actor.SetMapper( mapper )
        actor.GetProperty().SetColor( params[1], params[2], params[3] )
        actor.GetProperty().SetOpacity( params[4] )
        actor.GetProperty().BackfaceCullingOff()
        actor.PickableOff()

        return actor

    def add_volume_rep(self, name, data, style, params, cmd_string):
        at = self.cur_at()
        if not hasattr(at, "volume_reps"):
            at.volume_reps = {}

        if name not in at.volume_reps:
            # prepare transformation
            transform = vtk.vtkTransform()
            m = transform.GetMatrix()
            c = at.get_cell()
            # scale so 0--1 spans entire cell vector.  Note that data storage is opposite 
            # of lattice vector ordering
            c[0,:] /= data.shape[2]
            c[1,:] /= data.shape[1]
            c[2,:] /= data.shape[0]
            # transformation matrix is transpose of cell matrix
            for i0 in range(3):
                for i1 in range(3):
                    m.SetElement(i1,i0,c[i0,i1])

            # swap axes if needed
            # in data array, axes are 1 and 2 (swapaxes), but these correspond to lattice vectors 0 and 1 (i_of)
            if m.Determinant() < 0:
                data = np.ascontiguousarray(np.swapaxes(data, 1, 2))
                i_of = [1,0,2]
                for i0 in range(3):
                    for i1 in range(3):
                        m.SetElement(i1,i0,c[i_of[i0],i1])

            img = self.array_to_image(data)
            at.volume_reps[name] = (img, transform, [])
        else:
            img = at.volume_reps[name][0]
            transform = at.volume_reps[name][1]

        if style == "isosurface":
            at.volume_reps[name][2].append((self.make_isosurface(img, transform, params), cmd_string))
        else:
            raise ValueError("add_volume_rep got unsupported style '{}'".format(style))

    def prep_for_atoms_write(self, ats=None):
        if ats is None:
            ats = self.at_list

        for at in ats:
            # all _vtk_commands have already been incorporated in, should be part of state and no longer needed
            try:
                del at.info["_vtk_commands"]
            except KeyError:
                pass

            # save bonds
            try:
                bonds = at.bonds
            except AttributeError:
                bonds = None
            if bonds:
                at.bonds.write_to_atoms_arrays()

            # save volume commands
            if len(at.volume_reps) > 0:
                if "_vtk_commands" not in at.info:
                    at.info["_vtk_commands"] = ""
                for volume_rep in at.volume_reps.values():
                    for (_, cmd) in volume_rep[2]:
                        at.info["_vtk_commands"] += ' '.join(cmd) + " ;"
                if at.info["_vtk_commands"].endswith(";"):
                    at.info["_vtk_commands"] = at.info["_vtk_commands"][:-1]


    def prep_after_atoms_read(self, ats=None):
        if ats is None:
            ats = self.at_list

        for at in ats:
            if "_vtk_bonds" in at.arrays:
                at.bonds = DavTKBonds(at, self.settings)
                at.bonds.read_from_atoms_arrays()
