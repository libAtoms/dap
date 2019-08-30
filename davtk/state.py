import sys, ase.io, math
import numpy as np
import vtk
import ase.neighborlist
import re
from vtk.util.vtkImageImportFromArray import vtkImageImportFromArray
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

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

def get_atom_type_list(settings, at):
    atom_type_field = settings["atom_type_field"]

    if atom_type_field == "Z":
        atom_type_list = [str(Z) for Z in at.get_atomic_numbers()]
    elif atom_type_field == "species":
        atom_type_list = [sp for sp in at.get_chemical_symbols()]
    else:
        atom_type_list = [str(val) for val in at.arrays[atom_type_field]]

    return atom_type_list

def get_atom_prop(settings, atom_type, i=None, arrays=None):
    prop = vtk.vtkProperty()
    type_data = settings["atom_types"][atom_type]

    # only set color for non-colormap, otherwise just default to white (maybe should use some point on colormap?)
    if type_data["color"] is not None:
        prop.SetColor(type_data["color"])

    return prop

def get_atom_radius(settings, atom_type, i, at=None):
    if settings["atom_types"][atom_type]["radius_field"] is not None:
        radius_field = settings["atom_types"][atom_type]["radius_field"][0]
        factor = settings["atom_types"][atom_type]["radius_field"][1]
        if isinstance(i, list):
            atom_type_list = i
            r = np.median([at.arrays[radius_field][ii] for ii in range(len(at)) if atom_type_list[ii] == atom_type])
        else:
            r = at.arrays[radius_field][i]
        r *= factor
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

        atom_type_list = get_atom_type_list(self.settings, self.at)

        if in_cutoff is None or len(in_cutoff) == 0: # fully auto
            max_cutoff = max([none_zero(self.settings["atom_types"][atom_type_list[i]]["bonding_radius"]) for i in range(len(self.at))])
            u_cutoff_min = lambda i1, i2 : 0.0
            u_cutoff_max = lambda i1, i2 : 0.5 * ( self.settings["atom_types"][atom_type_list[i1]]["bonding_radius"] +
                                                   self.settings["atom_types"][atom_type_list[i2]]["bonding_radius"] )
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
                if d > 0.0 and d >= u_cutoff_min(i,j) and d <= u_cutoff_max(i, j) and pair_match(atom_type_list, i, j, at_type, at_type2):
                    self.bonds[i].append({ "j" : j, "S" : np.array(S), "name" : name, "picked" : False})
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
                m = re.search('^([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_(.*)$', bond_str)
                if m:
                    (j, S, picked, name) = ( int(m.group(1)), [int(m.group(i)) for i in range(2,5)], 
                                             str_to_bool(m.group(5)), m.group(6) )
                    if name not in self.settings["bond_types"]:
                        raise ValueError("Unknown bond_type '{}' reading from Atoms object, known types: {}\n".format(name, list(self.settings["bond_types"].keys()))+
                                         "            Perhaps you forgot to restore saved settings?")
                    self.bonds[at_i].append( {"j" : j, "S" : np.array(S), "name" : name, "picked" : picked } )

    def pair_mic(self, name, ind1, ind2):
        v = self.at.get_distance(ind1, ind2, mic=True, vector=True)
        self.bonds[ind1].append( {"j" : ind2, "S" : np.array([0,0,0]), "name" : name, "picked" : False } )
        self.bonds[ind2].append( {"j" : ind1, "S" : np.array([0,0,0]), "name" : name, "picked" : False } )

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
        self.shapes = self.create_shapes()
        self.create_vtk_structures()
        self.cur_frame = 0
        self.renderer = renderer
        self.iRen = iRen

        self.label_actor_pool = []
        self.cur_n_label_actors = 0

        self.bonds_points_data = vtk.vtkPoints()
        self.vectors_points_data = vtk.vtkPoints()
        self.atoms_points_data = {}

        self.atoms_actors = []
        self.image_atom_actors = []
        self.bonds_actor = None
        self.vector_actors = []
        self.volume_reps_actors = []

        self.legend_sphere_actors = []
        self.legend_label_actors = []

        self.cell_box_actor = None

        self.frame_label_actor = None

        self.saved_views = {}

        self.active = False

    def __len__(self):
        return len(self.at_list)

    def cur_at(self):
        return self.at_list[self.cur_frame]

    def create_shapes(self):
        sources = {}
        mappers = {}

        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1.0)
        sphere.SetPhiResolution(8)
        sphere.SetThetaResolution(16)
        sphere_mapper = vtk.vtkPolyDataMapper()
        sphere_mapper.SetInputConnection(sphere.GetOutputPort())
        sphere_mapper.Update()
        sources["sphere"] = (sphere, sphere_mapper.GetInput().GetNumberOfPoints())
        mappers["sphere"] = sphere_mapper

        cyl_y = vtk.vtkCylinderSource()
        cyl_y.SetRadius(1.0)
        cyl_y.SetHeight(1.0)
        cyl_y.SetResolution(8)
        cyl_y_mapper = vtk.vtkPolyDataMapper()
        cyl_y_mapper.SetInputConnection(cyl_y.GetOutputPort())
        cyl_y_mapper.Update()
        sources["cylinder_y"] = (cyl_y, cyl_y_mapper.GetInput().GetNumberOfPoints())
        mappers["cylinder_y"] = cyl_y_mapper

        arrow_x = vtk.vtkArrowSource()
        arrow_x.SetShaftRadius(1.0)
        arrow_x.SetTipRadius(2.5)
        arrow_x.SetTipLength(0.25)
        arrow_x.SetShaftResolution(8)
        arrow_x.SetTipResolution(16)
        arrow_x_mapper = vtk.vtkPolyDataMapper()
        arrow_x_mapper.SetInputConnection(arrow_x.GetOutputPort())
        arrow_x_mapper.Update()
        mappers["arrow_x"] = arrow_x_mapper

        cone_x = vtk.vtkConeSource()
        cone_x.SetRadius(1.0)
        cone_x.SetHeight(1.0)
        cone_x.SetResolution(16)
        cone_x.SetDirection(1.0, 0.0, 0.0)
        sources["cone_x"] = (cone_x,)

        cyl_x = vtk.vtkTransformPolyDataFilter()
        t = vtk.vtkTransform()
        t.RotateZ(-90)
        t.Translate(0.0, 0.5, 0.0)
        cyl_x.SetTransform(t)
        cyl_x.SetInputConnection(cyl_y.GetOutputPort())
        cyl_x_mapper = vtk.vtkPolyDataMapper()
        cyl_x_mapper.SetInputConnection(cyl_x.GetOutputPort())
        cyl_x_mapper.Update()
        mappers["cylinder_x_end_origin"] = cyl_x_mapper

        return { 'sources' : sources, 'mappers' : mappers }
        return sources

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

    def update_settings(self):
        self.create_bond_lut()
        self.create_atom_type_luts()

    def set_frame(self, frame):
        if frame.startswith("+") or frame.startswith("-"):
            self.cur_frame += int(frame)
        else:
            self.cur_frame = int(frame)
        self.cur_frame = self.cur_frame % len(self.at_list)

    def update(self, what=None):
        if not self.active:
            return

        if what is None or what == "cur":
            what = "+0"
        if not re.search('^(rotate|settings|[+-]?\d+)$', what):
            raise ValueError("update what='{}', not rotate or settings or number or [+|-]number".format(what))

        if what == "rotate":
            self.update_legend(self.cur_at())
            return

        if what != "settings":
            self.set_frame(what)
        at = self.cur_at()

        self.update_settings()

        self.renderer.SetBackground(self.settings["background_color"])

        self.update_cell_box(at.get_cell(), what == "settings")
        self.update_atoms(at, what == "settings")
        self.update_labels(at, what == "settings")
        self.update_bonds(at, what == "settings")
        self.update_vectors(at, what == "settings")
        self.update_image_atoms(at, what == "settings")
        self.update_volume_reps(at, what == "settings")
        self.update_legend(at, what == "settings")
        self.update_frame_label(what == "settings")

        # refresh display
        self.renderer.GetRenderWindow().Render()

    def update_volume_reps(self, at, settings_only = False):
        # clean out existing actors
        for actor in self.volume_reps_actors:
            self.renderer.RemoveActor(actor)
        self.volume_reps_actors = []

        # create new ones if needed
        if hasattr(at, "volume_reps"):
            for volume_rep in at.volume_reps.values():
                for (actor, _) in volume_rep[2]:
                    self.volume_reps_actors.append(actor)
                    self.renderer.AddActor(actor)

    def update_cell_box(self, cell, settings_only=False):
        if self.cell_box_actor is None:
            self.cell_box_actor = vtk.vtkActor()
            self.cell_box_actor._vtk_type = "cell_box"

        actor = self.cell_box_actor
        actor.GetProperty().SetColor(self.settings["cell_box_color"])
        actor.PickableOff()

        if settings_only:
            return

        pts = vtk.vtkPoints()
        for i0 in range(2):
            for i1 in range(2):
                for i2 in range(2):
                    pts.InsertNextPoint(i0*cell[0] + i1*cell[1] + i2*cell[2])
        # 0 0 0    0 0 1    0 1 0    0 1 1     1 0 0  1 0 1   1 1 0   1 1 1
        lines = vtk.vtkCellArray()
        segment_point_ids = [ ((0,0),(1,1)) , ((0,0),(1,2)) , ((0,0),(1,4)) , ((0,1),(1,3)) , ((0,1),(1,5)) , ((0,2),(1,3)) ,
                      ((0,2),(1,6)) , ((0,4),(1,5)) , ((0,4),(1,6)) , ((0,3),(1,7)) , ((0,5),(1,7)) , ((0,6),(1,7)) ]
        for segment in segment_point_ids:
            l = vtk.vtkLine()
            l.GetPointIds().SetId(segment[0][0], segment[0][1])
            l.GetPointIds().SetId(segment[1][0], segment[1][1])
            lines.InsertNextCell(l)

        linesPolyData = vtk.vtkPolyData()
        linesPolyData.SetPoints(pts)
        linesPolyData.SetLines(lines)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(linesPolyData)

        actor.SetMapper(mapper)
        self.renderer.AddActor(actor)

    def atoms_plotting_info(self, at):
        atom_type_list = get_atom_type_list(self.settings, at)
        pos = at.positions

        # create structures for each atom type
        points_lists = {}
        radius_lists = {}
        colormap_vals_lists = {}
        i_at_lists = {}
        unique_types = set([at_type for at_type in atom_type_list])
        for at_type in unique_types:
            points_lists[at_type] = []
            radius_lists[at_type] = []
            colormap_vals_lists[at_type] = []
            i_at_lists[at_type] = []

        # create position, radius, and colormap value lists, separately for each atom type
        for i_at in range(len(at)):
            at_type = atom_type_list[i_at]
            picked = at.arrays["_vtk_picked"][i_at]
            r = get_atom_radius(self.settings, at_type, i_at, at)

            colormap = self.settings["atom_types"][at_type]["colormap"]
            if colormap is None: # fixed color
                if picked:
                    colormap_val = (1.0, 0.0, 0.0)
                else:
                    colormap_val = (0.0, 0.0, 0.0)
            else:
                (colormap_name, colormap_field) = colormap
                if picked: # every colormap has picked_color at max range + 1.0
                    colormap_val = (self.settings["colormaps"][colormap_name][-1][0]+1.0, 0.0, 0.0)
                else:
                    colormap_val = (np.linalg.norm(at.arrays[colormap_field]), 0.0, 0.0)

            points_lists[at_type].append(pos[i_at])
            radius_lists[at_type].append(r)
            colormap_vals_lists[at_type].append(colormap_val)
            i_at_lists[at_type].append(i_at)

        return (unique_types, points_lists, radius_lists, colormap_vals_lists, i_at_lists)

    # need to see what can be optimized if settings_only is True
    def update_atoms(self, at, settings_only = False):

        # get structures 
        (unique_types, points_lists, radius_lists, colormap_vals_lists, i_at_lists) = self.atoms_plotting_info(at)

        # cleaup up actors
        for actor in self.atoms_actors:
            self.renderer.RemoveActor(actor)
        self.atoms_actors = []

        for at_type in unique_types:
            # create points
            if at_type not in self.atoms_points_data:
                self.atoms_points_data[at_type] = vtk.vtkPoints()
            points = self.atoms_points_data[at_type]
            points.Reset()
            points.SetNumberOfPoints(len(points_lists[at_type]))
            for pt_i in range(len(points_lists[at_type])):
                points.SetPoint(pt_i, points_lists[at_type][pt_i])

            # create polydata from points, radii, and values for colormap
            glyphs_data = vtk.vtkPolyData()
            glyphs_data.SetPoints(points)
            glyphs_data.GetPointData().SetScalars(numpy_to_vtk(np.array(radius_lists[at_type])))
            glyphs_data.GetPointData().SetVectors(numpy_to_vtk(np.array(colormap_vals_lists[at_type])))

            # create spherical glyphs
            glyphs = vtk.vtkGlyph3D()
            glyphs.SetInputData(glyphs_data)
            glyphs.SetSourceConnection(self.shapes["sources"]["sphere"][0].GetOutputPort())
            glyphs.SetScaleModeToScaleByScalar()
            glyphs.SetColorModeToColorByVector()
            glyphs.Update()

            id_glyphs = vtk.vtkIdFilter()
            id_glyphs.SetInputConnection(glyphs.GetOutputPort())
            id_glyphs.SetIdsArrayName("IDs")
            id_glyphs.FieldDataOn()

            # create mapper and set color LUT
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(id_glyphs.GetOutputPort())
            mapper.SetLookupTable(self.atom_type_luts[at_type])

            # create actor
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor._vtk_type = "atoms_glyphs"
            actor.point_to_input_point = self.shapes["sources"]["sphere"][1]
            actor.i_at = i_at_lists[at_type]

            # add and save actor
            self.renderer.AddActor(actor)
            self.atoms_actors.append(actor)

    # need to see what can be optimized if settings_only is True
    def update_vectors(self, at, settings_only = False):
        for actor in self.vector_actors:
            self.renderer.RemoveActor(actor)
        self.vector_actors = []

        if at.info.get("_vtk_vectors", None) is None:
            return

        pos = at.get_positions()

        rad = at.info["_vtk_vectors"]["radius"]
        vector_color = at.info["_vtk_vectors"]["color"]
        if vector_color == "atom":
            atom_type_list = get_atom_type_list(self.settings, at)

        if at.info["_vtk_vectors"]["field"] == "magmoms":
            vectors = at.get_magnetic_moments()
        elif at.info["_vtk_vectors"]["field"] == "initial_magmoms":
            vectors = at.get_initial_magnetic_moments()
        else:
            vectors = at.arrays[at.info["_vtk_vectors"]["field"]]

        # start out with fake orientation
        self.renderer.GetActiveCamera().OrthogonalizeViewUp()
        orientation = [ 0.0, 1.0, 0.0 ]

        # use vector field shape to determine whether orientation is meaningful
        if len(vectors.shape) == 2 and vectors.shape[1] == 3: # 3-vectors
            if vector_color == "sign":
                raise ValueError("Can't color vectors by sign when they are not actually scalars")
            vectors_use = vectors
            orientation = None
        elif len(vectors.shape) == 1: # scalars
            vectors_use = np.outer(vectors, orientation)
        else:
            raise ValueError("Don't know how to draw vectors for field with shape {}".format(vectors.shape))

        vectors_use *= at.info["_vtk_vectors"]["scale"]

        vector_norms = np.linalg.norm(vectors_use, axis=1)
        vectors_hat = (vectors_use.T / vector_norms).T

        cyl_x_mapper = self.shapes["mappers"]["cylinder_x_end_origin"]
        arrow_x_mapper = self.shapes["mappers"]["arrow_x"]
        cone_x_source = self.shapes["sources"]["cone_x"][0]
        source_orient = [1.0, 0.0, 0.0]

        for i_at in range(len(at)):

            axis = np.cross(vectors_hat[i_at], source_orient)
            if np.linalg.norm(axis) < 1.0e-4: # essentially parallel to source_orient,  set angle to 0, axis irrelevant
                axis = (1.0,0.0,0.0)
                angle = 0.0
            else:
                axis /= np.linalg.norm(axis)
                angle = -np.arccos(np.dot(vectors_hat[i_at], source_orient))*180.0/np.pi

            if vector_color == "atom":
                at_type = atom_type_list[i_at] 
                color = get_atom_prop(self.settings, at_type, i_at, at.arrays).GetColor()
            elif vector_color == "sign":
                if vectors[i_at] < 0:
                    color = at.info["_vtk_vectors"]["sign_colors"][3:6]
                else:
                    color = at.info["_vtk_vectors"]["sign_colors"][0:3]
            else:
                color = at.info["_vtk_vectors"]["color"]

            if orientation is not None:
                actors = [vtk.vtkFollower(),vtk.vtkFollower(),vtk.vtkFollower()]
                for actor in actors:
                    actor.SetCamera(self.renderer.GetActiveCamera())
                    actor.SetPosition(pos[i_at])
                    # assumes source_orient = \hat{x}

                # 0 gets arrowhead
                actors[0].SetMapper(cyl_x_mapper)
                actors[0].RotateWXYZ(angle, axis[0], axis[1], axis[2])
                actors[0].SetScale(vector_norms[i_at]/2.0-3.0*rad, rad, rad)
                actors[1].SetMapper(cyl_x_mapper)
                actors[1].RotateWXYZ(angle+180.0, axis[0], axis[1], axis[2])
                actors[1].SetScale(vector_norms[i_at]/2.0, rad, rad)

                cone = vtk.vtkTransformPolyDataFilter()
                t = vtk.vtkTransform()
                t.Translate((vector_norms[i_at]/2.0-1.5*rad) / (3.0*rad), 0.0, 0.0)
                cone.SetTransform(t)
                cone.SetInputConnection(cone_x_source.GetOutputPort())
                cone_mapper = vtk.vtkPolyDataMapper()
                cone_mapper.SetInputConnection(cone.GetOutputPort())
                actors[2].SetMapper(cone_mapper)
                actors[2].RotateWXYZ(angle, axis[0], axis[1], axis[2])
                actors[2].SetScale(rad*3.0, rad*2.5, rad*2.5)

            else:
                actors = [vtk.vtkActor()]
                actors[0].SetMapper(arrow_x_mapper)
                # assumes source_orient = \hat{x}
                actors[0].SetScale(vector_norms[i_at], rad, rad)
                actors[0].SetPosition(pos[i_at]-vector_norms[i_at]/2.0*vectors_hat[i_at])
                actors[0].RotateWXYZ(angle, axis[0], axis[1], axis[2])

            for actor in actors:
                actor._vtk_type = "vector"
                actor.VisibilityOn()
                actor.PickableOff()
                actor.GetProperty().SetColor(color)
                self.renderer.AddActor(actor)
                self.vector_actors.append(actor)

    # need to see what can be optimized if settings_only is True
    def update_bonds(self, at, settings_only = False):
        if not hasattr(at, "bonds"):
            if self.bonds_actor is not None:
                self.bonds_actor.VisibilityOff()
                self.bonds_actor.PickableOff()
                self.renderer.RemoveActor(self.bonds_actor)
            return

        if self.bonds_actor is None:
            self.bonds_actor = vtk.vtkActor()
            self.bonds_actor._vtk_type = "bonds_glyphs"
            self.bonds_actor.point_to_input_point = self.shapes["sources"]["cylinder_y"][1]

        actor = self.bonds_actor

        pos = at.get_positions()

        i_at_bond = []

        points_list = []
        axes_list = []
        angles_list = []
        scales_list = []
        color_indices_list = []

        cell = at.get_cell()
        for i_at in range(len(at)):
            for (i_bond, b) in enumerate(at.bonds[i_at]):
                j_at = b["j"]
                if i_at < j_at:
                    continue
                dr = pos[j_at] - pos[i_at] + np.dot(b["S"], cell)
                dr_norm = np.linalg.norm(dr)
                picked = b["picked"]
                name = b["name"]
                rad = self.settings["bond_types"][name]["radius"]

                axis = np.cross(dr, [0.0, 1.0, 0.0])
                if np.linalg.norm(axis) < 1.0e-4: # essentially parallel to \hat{y}
                    axis = (1.0,0.0,0.0)
                    angle = 0.0
                else:
                    axis /= np.linalg.norm(axis)
                    dr_hat = dr / dr_norm
                    # angle = -np.arccos(np.dot(dr_hat, [0.0, 1.0, 0.0]))*180.0/np.pi
                    angle = -np.arccos(dr_hat[1])*180.0/np.pi

                color_index = 0 if picked else self.settings["bond_types"][name]["index"]

                points_list.append(pos[i_at]+dr/4.0)
                axes_list.append(axis)
                angles_list.append(angle)
                scales_list.append((rad, dr_norm/2.0, rad))
                color_indices_list.append((color_index,))
                i_at_bond.append((i_at, i_bond))

                points_list.append(pos[j_at]-dr/4.0)
                axes_list.append(axis)
                angles_list.append(angle)
                scales_list.append((rad, dr_norm/2.0, rad))
                color_indices_list.append((color_index,))
                i_at_bond.append((i_at, i_bond))

        if len(points_list) == 0:
            self.bonds_actor.VisibilityOff()
            self.bonds_actor.PickableOff()
            self.renderer.RemoveActor(self.bonds_actor)
            return

        # data needed to place, rotate, and scale glyphs
        points = self.bonds_points_data
        points.Reset()
        points.SetNumberOfPoints(len(points_list))
        for pt_i in range(len(points_list)):
            points.SetPoint(pt_i, points_list[pt_i])
        axes = numpy_to_vtk(np.array(axes_list))
        angles = numpy_to_vtk(np.array(angles_list))
        scales = numpy_to_vtk(np.array(scales_list))
        color_indices = numpy_to_vtk(np.array(color_indices_list, dtype=int))
        color_indices.SetName("color_indices")

        glyphs_data = vtk.vtkPolyData()
        glyphs_data.SetPoints(points)
        glyphs_data.GetPointData().AddArray(axes)
        glyphs_data.GetPointData().AddArray(angles)
        glyphs_data.GetPointData().AddArray(scales)
        glyphs_data.GetPointData().SetScalars(color_indices)

        def place_glyph(__vtk__temp0=0,__vtk__temp1=0):
            point = glyphs.GetPoint()
            point_id = glyphs.GetPointId()
            d = glyphs.GetPointData()
            axis =  d.GetAbstractArray(0).GetTuple(point_id)
            angle = d.GetAbstractArray(1).GetTuple(point_id)
            scale = d.GetAbstractArray(2).GetTuple(point_id)

            trans = vtk.vtkTransform()
            trans.Translate(point)
            trans.RotateWXYZ(angle[0], axis)
            trans.Scale(scale[0], scale[1], scale[2])

            trans_cyl = vtk.vtkTransformPolyDataFilter()
            trans_cyl.SetInputConnection(self.shapes["sources"]["cylinder_y"][0].GetOutputPort())
            trans_cyl.SetTransform(trans)
            glyphs.SetSourceConnection(trans_cyl.GetOutputPort())

        glyphs = vtk.vtkProgrammableGlyphFilter()
        glyphs.SetSourceConnection(self.shapes["sources"]["cylinder_y"][0].GetOutputPort())
        glyphs.SetInputData(glyphs_data)
        glyphs.SetGlyphMethod(place_glyph)
        glyphs.SetColorModeToColorByInput()

        id_glyphs = vtk.vtkIdFilter()
        id_glyphs.SetInputConnection(glyphs.GetOutputPort())
        id_glyphs.SetIdsArrayName("IDs")
        id_glyphs.FieldDataOn() # check Off?

        glyphs_mapper = vtk.vtkPolyDataMapper()
        glyphs_mapper.SetInputConnection(id_glyphs.GetOutputPort())

        glyphs_mapper.SetLookupTable(self.bond_lut)

        actor.SetMapper(glyphs_mapper)
        actor.i_at_bond = i_at_bond

        actor.VisibilityOn()
        actor.PickableOn()
        self.renderer.AddActor(actor)

    def update_image_atoms(self, at, settings_only = False):
        if "_vtk_images" not in at.info:
            return

        (unique_types, points_lists, radius_lists, colormap_vals_lists, i_at_lists) = self.atoms_plotting_info(at)

        n_images = 0
        pos = at.get_positions()
        cell = at.get_cell()
        cell_inv = at.get_reciprocal_cell().T
        for at_type in unique_types:
            for (point, r, colormap_val, i_at) in zip(points_lists[at_type], radius_lists[at_type], colormap_vals_lists[at_type], i_at_lists[at_type]):
                n0 = int(math.ceil(at.info["_vtk_images"][0]))
                n1 = int(math.ceil(at.info["_vtk_images"][1]))
                n2 = int(math.ceil(at.info["_vtk_images"][2]))
                for i0 in range(-n0, n0+1):
                    for i1 in range(-n1, n1+1):
                        for i2 in range(-n2, n2+1):
                            if (i0,i1,i2) == (0,0,0):
                                continue
                            image_pos = pos[i_at] + np.dot([i0, i1, i2], cell)
                            image_pos_scaled = np.dot(image_pos, cell_inv)
                            if (image_pos_scaled[0] >= -at.info["_vtk_images"][0] and
                                image_pos_scaled[0] <= 1+at.info["_vtk_images"][0] and
                                image_pos_scaled[1] >= -at.info["_vtk_images"][1] and
                                image_pos_scaled[1] <= 1+at.info["_vtk_images"][1] and
                                image_pos_scaled[2] >= -at.info["_vtk_images"][2] and
                                image_pos_scaled[2] <= 1+at.info["_vtk_images"][2]):
                                n_images += 1
                                if n_images > len(self.image_atom_actors):
                                    image_actor = vtk.vtkActor()
                                    self.image_atom_actors.append(image_actor)
                                    self.renderer.AddActor(image_actor)
                                    image_actor.VisibilityOn()
                                    image_actor.PickableOn()
                                    image_actor._vtk_type = "image_atom"
                                    image_actor.SetMapper(self.shapes["mappers"]["sphere"])
                                else:
                                    image_actor = self.image_atom_actors[n_images-1]

                                color = self.atom_type_luts[at_type].GetColor(np.linalg.norm(colormap_val))
                                image_actor.GetProperty().SetColor(color)
                                image_actor.SetPosition(image_pos)
                                image_actor.SetScale(r,r,r)
                                image_actor.i_at = i_at
                                image_actor.VisibilityOn()

        if len(self.image_atom_actors) > n_images:
            for actor in self.image_atom_actors[n_images:]:
                actor.VisibilityOff()
                actor.PickableOff()
                actor.i_at = -1

    def update_legend(self, at, settings_only = False):

        # remove, free, and re-create (recycling seems to suffer from bug perhaps addressed
        # in merge request 3904)
        for actor in self.legend_sphere_actors + self.legend_label_actors:
            self.renderer.RemoveActor(actor)
        self.legend_sphere_actors = []
        self.legend_label_actors = []

        if not self.settings["legend"]["show"]:
            return

        pos = at.get_positions()

        display_size = self.renderer.GetRenderWindow().GetSize()
        dp_world = self.label_offset_world(pos[0])

        atom_type_list = get_atom_type_list(self.settings, at)
        unique_atom_types = sorted(list(set(atom_type_list)))

        for i in range(len(unique_atom_types)):
            self.legend_sphere_actors.append(vtk.vtkActor())
            self.legend_sphere_actors[-1]._vtk_type = "legend_sphere"
            self.legend_label_actors.append(vtk.vtkBillboardTextActor3D())
            self.legend_label_actors[-1]._vtk_type = "legend_label"

        max_r = 0.0
        longest_label = ""
        for (l_i, (legend_sphere_actor, legend_label_actor, at_type)) in enumerate(zip(self.legend_sphere_actors, self.legend_label_actors, unique_atom_types)):
            max_r = max(max_r, get_atom_radius(self.settings, at_type, list(atom_type_list), at))
            longest_label = longest_label if len(longest_label) > len(at_type) else at_type
        legend_raw_pos = self.settings["legend"]['position']
        spacing = self.settings["legend"]["spacing"] * max(int(2.5*max_r/dp_world), self.settings["atom_label"]["fontsize"])
        legend_offset = np.zeros((2),dtype=int)
        if legend_raw_pos[0] >= 0:
            legend_offset[0] = int(1.3 * max_r/dp_world)
        else:
            legend_offset[0] = -int( 1.3 * max_r/dp_world + len(longest_label)*self.settings["atom_label"]["fontsize"]/1.5 )
        if legend_raw_pos[1] >= 0:
            legend_offset[1] = int(spacing * (len(unique_atom_types)-0.5))
        else:
            legend_offset[1] = -int(0.5*spacing)
        legend_pos = (legend_raw_pos + legend_offset) % display_size

        for (l_i, (legend_sphere_actor, legend_label_actor, at_type)) in enumerate(zip(self.legend_sphere_actors, self.legend_label_actors, unique_atom_types)):
            atom_i_of_type = atom_type_list.index(at_type)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(self.shapes["sources"]["sphere"][0].GetOutputPort())
            legend_sphere_actor.SetMapper(mapper)

            prop = get_atom_prop(self.settings, at_type, atom_i_of_type, at.arrays)

            r = get_atom_radius(self.settings, at_type, list(atom_type_list), at)
            legend_sphere_actor.SetProperty(prop)
            legend_sphere_actor.SetScale(r, r, r)

            self.renderer.SetWorldPoint(list(pos[0]) + [1.0])
            self.renderer.WorldToDisplay()
            arb_point = self.renderer.GetDisplayPoint()

            self.renderer.SetDisplayPoint(legend_pos[0], legend_pos[1]-spacing*l_i, arb_point[2])
            self.renderer.DisplayToWorld()
            sphere_pos_world = self.renderer.GetWorldPoint()
            sphere_pos_world = np.array(sphere_pos_world[0:3]) / sphere_pos_world[3]

            legend_sphere_actor.SetPosition(sphere_pos_world)
            legend_sphere_actor.PickableOff()
            legend_sphere_actor.VisibilityOn()

            legend_label_actor.SetInput(at_type)
            legend_label_actor.SetPosition(sphere_pos_world)
            if dp_world > 0:
                max_r_disp = max_r/dp_world
                text_disp = self.settings["atom_label"]["fontsize"]
            else:
                max_r_disp = 0
            legend_label_actor.SetDisplayOffset(int(1.3*max_r_disp), -int(text_disp/2.0))
            legend_label_actor.SetTextProperty(self.settings["atom_label"]["prop"])
            legend_label_actor.PickableOff()
            legend_label_actor.VisibilityOn()

        for actor in self.legend_sphere_actors + self.legend_label_actors:
            self.renderer.AddActor(actor)

    def update_frame_label(self, settings_only = False):
        if self.frame_label_actor is None:
            self.frame_label_actor = vtk.vtkTextActor()
            self.frame_label_actor._vtk_type = "frame_label"
            self.frame_label_actor.PickableOff()

        self.frame_label_actor.SetTextProperty(self.settings["frame_label"]["prop"])
        self.frame_label_actor.SetDisplayPosition(20,20)

        if self.settings["frame_label"]["field"] == "_NONE_":
            frame_label_str = ""
        elif self.settings["frame_label"]["field"] == "config_n":
            frame_label_str = str(self.cur_frame)
        else:
            frame_label_str = str(self.cur_at().info.get(self.settings["frame_label"]["field"],"None"))

        self.frame_label_actor.SetInput(frame_label_str)
        self.renderer.AddActor(self.frame_label_actor)

    def vector_lut(self, colors):
        # look up table for bond colors
        vector_lut = vtk.vtkColorTransferFunction()
        vector_lut.SetRange(0, len(colors)-1)
        for (i, c) in enumerate(colors):
            vector_lut.AddRGBPoint(float(i), c[0], c[1], c[2])
        return vector_lut

    def create_bond_lut(self):
        # look up table for bond colors
        self.bond_lut = vtk.vtkColorTransferFunction()
        n_types = len(self.settings["bond_types"])
        self.bond_lut.SetRange(0, n_types)
        self.bond_lut.AddRGBPoint(0.0, self.settings["picked_color"][0], self.settings["picked_color"][1], self.settings["picked_color"][2])
        for i in range(1, n_types+1):
            name = self.settings["bond_name_of_index"][i]
            self.bond_lut.AddRGBPoint(float(i), self.settings["bond_types"][name]["color"][0],
                                      self.settings["bond_types"][name]["color"][1],
                                      self.settings["bond_types"][name]["color"][2])

    def create_atom_type_luts(self):
        # make sure every atom type exists, autogenerating if needed
        for at_type in sorted(list(set(get_atom_type_list(self.settings, self.cur_at())))):
            s = self.settings["atom_types"][at_type]

        self.luts = {}
        # create true color maps
        for (name, data) in self.settings["colormaps"].items():
            lut = vtk.vtkColorTransferFunction()
            lut.SetRange(data[0][0], data[-1][0])
            for pt in data:
                lut.AddRGBPoint(pt[0], pt[1], pt[2], pt[3])
            lut.AddRGBPoint(data[-1][0]+1.0, self.settings["picked_color"][0], self.settings["picked_color"][1], self.settings["picked_color"][2])
            self.luts[name] = lut

        self.atom_type_luts = {}
        # create trivial color maps for fixed colors, or link to true maps
        for name in self.settings["atom_types"].get_all():
            color = self.settings["atom_types"][name]["color"]
            if color is not None: # just a color
                lut = vtk.vtkColorTransferFunction()
                lut.SetRange(0,0)
                lut.AddRGBPoint(0.0, color[0], color[1], color[2])
                lut.AddRGBPoint(1.0, self.settings["picked_color"][0], self.settings["picked_color"][1], self.settings["picked_color"][2])
                self.atom_type_luts[name] = lut
            elif self.settings["atom_types"][name]["colormap"] is not None:
                self.atom_type_luts[name] = self.luts[self.settings["atom_types"][name]["colormap"][0]]
            else:
                self.atom_type_luts[name] = None

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

    # need to see what can be optimized if settings_only is True
    # can/should this be done with some sort of GlyphMapper?
    def update_labels(self, at, settings_only = False):
        if not at.info["_vtk_show_labels"]:
            for actor in self.label_actor_pool:
                actor.SetVisibility(False)
            return

        if len(at) > len(self.label_actor_pool):
            prev_pool_size = len(self.label_actor_pool)
            new_actors = [vtk.vtkBillboardTextActor3D() for i in range(len(at)-prev_pool_size)]
            for actor in new_actors:
                actor._vtk_type="atom_label"
                self.renderer.AddActor(actor)
                actor.PickableOff()
            self.label_actor_pool.extend(new_actors)

        atom_type_list = get_atom_type_list(self.settings, at)
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
                label_field = self.settings["atom_label"]["field"]
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
            r = get_atom_radius(self.settings, atom_type_list[i_at], i_at, at)
            if dp_world > 0:
                dp_disp = 0.7*r/dp_world
            else:
                dp_disp = 0
            label_actor.SetDisplayOffset(int(dp_disp), int(dp_disp))
            label_actor.SetTextProperty(self.settings["atom_label"]["prop"])
            label_actor.SetVisibility(True)

        for actor in self.label_actor_pool[len(at):]:
            actor.SetVisibility(False)
        self.cur_n_label_actors = len(at)

    def create_vtk_structures(self, frames=None):
        for frame_i in self.frame_list(frames):
            at = self.at_list[frame_i]

            at.arrays["_vtk_picked"] = np.array([False] * len(at))
            at.info["_vtk_show_labels"] = False

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
        self.active = True
        self.update(what='0')

    def array_to_image(self, data):
        img = vtkImageImportFromArray()
        img.SetArray(data)
        img.Update()
        return img

    def delete_volume_rep(self, name):
        at = self.cur_at()
        for actor in at.volume_reps[name][2]:
            self.renderer.RemoveActor(actor)
        del at.volume_reps[name]

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
        actor._vtk_type = "isosurface"
        actor.SetMapper( mapper )
        actor.GetProperty().SetColor( params[1], params[2], params[3] )
        actor.GetProperty().SetOpacity( params[4] )
        actor.GetProperty().BackfaceCullingOff()
        actor.PickableOff()

        return actor

    def add_volume_rep(self, name, data, style, params, cmd_string):
        # at.volume_reps[name] is 3-element tuple: image, transform (only one for each data set), 
        #     and list of actors (could be multiple views of same data).
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
