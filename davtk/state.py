import sys, ase.io, math
import numpy as np
import vtk
import ase.neighborlist

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
    if "_vtk_type" in at.arrays:
        atom_type = at.arrays["_vtk_type"]
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
        r = arrays[settings["atom_types"][atom_type]["radius_field"]][i]
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

        def pair_match(at_type_a, i, j, at_type, at_type2):
            return ( ((at_type == '*' or str(at_type_a[i]) == at_type) and (at_type2 == '*' or str(at_type_a[j]) == at_type2)) or
                     ((at_type == '*' or str(at_type_a[j]) == at_type) and (at_type2 == '*' or str(at_type_a[i]) == at_type2)) )

        atom_type_array = get_atom_type_a(self.at)

        if in_cutoff is None or len(in_cutoff) == 0: # fully auto
            max_cutoff = max([self.settings["atom_types"][atom_type_array[i]]["bonding_radius"] for i in range(len(self.at))])
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
            if d > 0.0 and d >= u_cutoff_min(i,j) and d <= u_cutoff_max(i, j) and pair_match(atom_type_array, i, j, at_type, at_type2):
                self.bonds[i].append({ "j" : j, "v" : v, "d" : d, "S" : S, "name" : name, "picked" : False})

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

        self.update()

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

    def update_labels(self, at):
        if len(at) > len(self.label_actor_pool):
            prev_pool_size = len(self.label_actor_pool)
            self.label_actor_pool.extend([vtk.vtkBillboardTextActor3D() for i in range(len(at)-prev_pool_size)])

        atom_type_array = get_atom_type_a(at)
        pos = at.get_positions()

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
                else:
                    label_str = str(at.arrays[label_field][i_at])
            label_actor.SetInput(label_str)
            label_actor.SetPosition(pos[i_at])
            r = get_atom_radius(self.settings, atom_type_array[i_at], i_at, at.arrays)
            if i_at == 0: # figure out mapping from screen to world distances
                self.renderer.GetActiveCamera()
                # screen pos of arbitrary world point
                pt_disp = [0.0, 0.0, 0.0]
                self.iRen.GetInteractorStyle().ComputeWorldToDisplay(self.renderer, pos[i_at][0], pos[i_at][1], pos[i_at][2], pt_disp)
                # world pos of point 10 pixel away
                pt_disp[0] += 10.0
                pt_world = [0.0, 0.0, 0.0, 0.0]
                self.iRen.GetInteractorStyle().ComputeDisplayToWorld(self.renderer, pt_disp[0], pt_disp[1], pt_disp[2], pt_world)
                dp_world = np.linalg.norm(np.array(pos[i_at])-np.array(pt_world[0:3]))/10
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

        atom_type_array = get_atom_type_a(at)
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

            actor = self.at_list[frame_i].info["_NOPRINT_vtk_cell_box_actor"]
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(self.settings["cell_box_color"])
            actor.PickableOff()

    def create_vtk_structures(self, frames=None):
        for frame_i in self.frame_list(frames):
            at = self.at_list[frame_i]

            at.arrays["_vtk_picked"] = np.array([False] * len(at))
            at.info["_NOPRINT_vtk_cell_box_actor"] = vtk.vtkActor()
            at.info["_vtk_show_labels"] = False

    def show_frame(self, dframe=None, frame_i=None):
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

        at = self.at_list[self.cur_frame]

        self.update_atoms(at)
        self.update_labels(at)
        self.update_bonds(at)

        # remove all existing actors
        self.renderer.RemoveAllViewProps()

        # actor for cell box
        self.renderer.AddActor(at.info["_NOPRINT_vtk_cell_box_actor"])

        # create actors for atoms
        pos = at.get_positions()
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

        # need to do other actors, e.g. labels and bonds
        for actor in self.bond_actor_pool[0:self.cur_n_bond_actors]:
            self.renderer.AddActor(actor)

        if "_vtk_show_labels" in at.info and at.info["_vtk_show_labels"]:
            for actor in self.label_actor_pool[0:self.cur_n_label_actors]:
                self.renderer.AddActor(actor)

        # config_n
        txt = vtk.vtkTextActor()
        txt.SetInput(str(self.cur_frame))
        txt.SetTextProperty(self.settings["config_n_text_prop"])
        txt.SetDisplayPosition(20,20)
        self.renderer.AddActor(txt)

        # refresh display
        self.renderer.GetRenderWindow().Render()

    def measure_picked(self):
        at = self.at_list[self.cur_frame]

        at_indices = np.where(at.arrays["_vtk_picked"])[0]
        p = at.get_positions()

        print("measure:")
        for i_at in at_indices:
            print("atom {} pos {} {} {}".format(i_at,p[i_at][0],p[i_at][1],p[i_at][2]))
        for i in range(len(at_indices)):
            i_at = at_indices[i]
            for j_at in at_indices[i+1:]:
                dv = -at.get_distance(i_at,j_at,mic=True,vector=True)
                print("atom-distance {} {} vec {} {} {} ({})".format(i_at,j_at,dv[0],dv[1],dv[2],np.linalg.norm(dv)))

        if hasattr(at, "bonds"):
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

    def snapshot(self, filename, mag=1):
        renderLarge = vtk.vtkRenderLargeImage()
        renderLarge.SetInput(self.renderer)
        renderLarge.SetMagnification(mag)

        writer = vtk.vtkPNGWriter()
        writer.SetInputConnection(renderLarge.GetOutputPort())
        writer.SetFileName(filename)
        writer.Write()
