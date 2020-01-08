import sys, ase.io, math
from warnings import warn
import numpy as np, scipy
import vtk
import ase.neighborlist
import re
from vtk.util.vtkImageImportFromArray import vtkImageImportFromArray
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from ase.build.supercells import make_supercell

def bond_vector(cell, pos, i_at, j_at, S, dist=True):
    D = pos[j_at] - pos[i_at] + np.dot(S, cell)
    if dist:
        return (D, np.linalg.norm(D))
    else:
        return D

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
        if radius_field not in at.arrays:
            raise RuntimeError("Atom radius field '{}' is not available".format(radius_field))
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

    def cutoff(self, name, in_cutoff, at_type, at_type2, across_pbc=False):

        def none_zero(x):
            return x if x is not None else 0.0

        def pair_type_match(at_type_a, i, j, at_type, at_type2):
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

        nn_list = ase.neighborlist.neighbor_list('ijdS', self.at, max_cutoff, self_interaction=True)
        for (i, j, d, S) in zip(nn_list[0], nn_list[1], nn_list[2], nn_list[3]):
            try:
                if d > 0.0 and d >= u_cutoff_min(i,j) and d <= u_cutoff_max(i, j) and pair_type_match(atom_type_list, i, j, at_type, at_type2) and (across_pbc or np.all(S == 0)):
                    if i == j and S in (b["S"] for b in self.bonds[i]): # don't create opposite bonds for i-i (periodic image)
                        continue
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
                v = [ b['j'], b['S'][0], b['S'][1], b['S'][2], b['picked'], b['name'] ]
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
                    self.bonds[at_i].append( {"j" : j, "S" : np.array(S), "name" : name, "picked" : picked } )

    def pair_mic(self, name, ind1, ind2):
        if ind1 == ind2:
            raise ValueError("Can't minimum-image-convention bond atom to itself")
        v = self.at.get_distance(ind1, ind2, mic=True, vector=True)
        self.bonds[ind1].append( {"j" : ind2, "S" : np.array([0,0,0]), "name" : name, "picked" : False } )
        self.bonds[ind2].append( {"j" : ind1, "S" : np.array([0,0,0]), "name" : name, "picked" : False } )

    def set_picked(self, i_at, j_ind, stat):
        b = self.bonds[i_at][j_ind]
        b["picked"] = stat
        j_at = b["j"]
        if j_at != i_at: # find opposite
            for (bb_i, bb) in enumerate(self.bonds[j_at]):
                if bb["j"] == i_at and np.all(b["S"] == -bb["S"]):
                    self.bonds[b["j"]][bb_i]["picked"] = stat
                    return
            raise ValueError("set_picked failed to find opposite for {} {}".format(i_at, j_at))

    def delete_one(self, i_at, j_ind):
        b = self.bonds[i_at][j_ind]
        del self.bonds[i_at][j_ind]
        j_at = b["j"]
        if i_at != j_at: # remove opposite
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
        self.cur_frame = 0
        self.renderer = renderer
        self.iRen = iRen

        self.atom_label_actors = []

        self.bond_prop = {}
        self.bonds_points_data = {}
        self.vectors_points_data = vtk.vtkPoints()
        self.atoms_points_data = {}

        self.atoms_actors = []
        self.image_atom_actors = []
        self.bonds_actors = {}
        self.vector_actors = []
        self.volume_reps_actors = []
        self.volume_rep_prop = {}

        self.legend_sphere_actors = []
        self.legend_label_actors = []

        self.polyhedra_prop = {}
        self.polyhedra_actors = []

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
        sphere.SetPhiResolution(16)
        sphere.SetThetaResolution(32)
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

    def delete(self, at, atom_selection=None, bond_selection=None):
        # delete atoms and their bonds
        if "_vtk_orig_indices" not in at.arrays:
            at.new_array("_vtk_orig_indices", np.array(range(len(at))))
        if atom_selection is not None:
            if atom_selection == "picked":
                if "_vtk_picked" in at.arrays:
                    at_inds = np.where(at.arrays["_vtk_picked"])[0]
                else:
                    at_inds = []
            elif isinstance(atom_selection,int):
                at_inds = np.array([atom_selection])
            else:
                at_inds = np.array(atom_selection)
            del at[at_inds]
            if hasattr(at, "bonds"):
                at.bonds.delete_atoms(at_inds)

        # delete requested bonds
        if bond_selection is not None and hasattr(at, "bonds"):
            if bond_selection == "picked":
                for i_at in range(len(at)):
                    while True:
                        b_row = at.bonds[i_at]
                        try:
                            j_ind = next(jj for jj in range(len(b_row)) if b_row[jj]["picked"])
                        except StopIteration: # no more matching neighbors
                            break
                        at.bonds.delete_one(i_at, j_ind)
            elif bond_selection == "all":
                at.bonds.reinit()
            else:
                for i_at in range(len(at)):
                    while True:
                        b_row = at.bonds[i_at]
                        try:
                            j_ind = next(jj for jj in range(len(b_row)) if b_row[jj]["name"] == bond_selection)
                        except StopIteration: # no more matching neighbors
                            break
                        at.bonds.delete_one(i_at, j_ind)

    def update_settings(self):
        self.create_atom_type_luts()

    def set_frame(self, frame):
        if frame.startswith("+") or frame.startswith("-"):
            self.cur_frame += int(frame)
        else:
            self.cur_frame = int(frame)
        self.cur_frame = self.cur_frame % len(self.at_list)

    def update(self, what=None):

        if what is None or what == "cur":
            what = "+0"
        if not re.search('^(rotate|settings|color_only|[+-]?\d+)$', what):
            raise ValueError("update what='{}', not 'rotate' or 'settings' or 'color_only' or int or [+|-]int".format(what))

        if what not in ["settings", "color_only", "rotate"]:
            self.set_frame(what)

        if not self.active:
            return

        if what == "rotate":
            self.update_legend(self.cur_at())
            self.update_atom_labels(self.cur_at())
            return

        at = self.cur_at()

        self.update_settings()

        self.renderer.SetBackground(self.settings["background_color"])

        self.update_cell_box(at.get_cell(), what == "settings" or what == "color_only")
        self.update_atom_spheres(at)
        self.update_atom_labels(at)
        self.update_bonds(at, what == "color_only")
        self.update_vectors(at)
        self.update_volume_reps(at, what == "settings" or what == "color_only")
        self.update_legend(at)
        self.update_frame_label(at)
        self.update_polyhedra(at, what == "settings" or what == "color_only")

        # refresh display
        self.renderer.GetRenderWindow().Render()

    def update_volume_reps(self, at, settings_only = False):
        if settings_only:
            return

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
            actor = vtk.vtkActor()
            actor._vtk_type = "cell_box"
            actor.PickableOff()
            actor.SetProperty(self.settings["cell_box"]["prop"])
            self.cell_box_actor = actor

        actor = self.cell_box_actor

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

    def visible_images(self, at):
        pos = at.positions
        cell = at.get_cell()
        cell_inv = at.get_reciprocal_cell().T

        vis_images = []
        if "_vtk_images" in at.info:
            rv = at.info["_vtk_images"]
            for i_at in range(len(at)):
                p_list = []
                vis_images.append(p_list)

                for s0 in range(int(math.floor(rv[0])), int(math.ceil(rv[1]))):
                    for s1 in range(int(math.floor(rv[2])), int(math.ceil(rv[3]))):
                        for s2 in range(int(math.floor(rv[4])), int(math.ceil(rv[5]))):

                            s = (s0, s1, s2)
                            p = pos[i_at] + np.dot(s, cell)
                            p_scaled = np.dot(p, cell_inv)
                            if (p_scaled[0] >= rv[0] and p_scaled[0] < rv[1] and
                                p_scaled[1] >= rv[2] and p_scaled[1] < rv[3] and
                                p_scaled[2] >= rv[4] and p_scaled[2] < rv[5]):
                                p_list.append((p, s))
        else:
            for i_at in range(len(at)):
                vis_images.append([(pos[i_at], (0,0,0))])

        return vis_images

    def atoms_plotting_info(self, at):
        atom_type_list = get_atom_type_list(self.settings, at)
        pos = at.positions
        cell = at.get_cell()
        cell_inv = at.get_reciprocal_cell().T
        rv = at.info.get("_vtk_images", (0.0, 1.0, 0.0, 1.0, 0.0, 1.0))

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

        vis_images = self.visible_images(at)

        # create position, radius, and colormap value lists, separately for each atom type
        for i_at in range(len(at)):
            at_type = atom_type_list[i_at]
            try:
                picked = at.arrays["_vtk_picked"][i_at]
            except KeyError:
                picked = False
            r = get_atom_radius(self.settings, at_type, i_at, at)

            colormap = self.settings["atom_types"][at_type]["colormap"]
            if colormap is None: # fixed color
                if picked:
                    colormap_val = (1.0, 0.0, 0.0)
                else:
                    colormap_val = (0.0, 0.0, 0.0)
            else:
                (colormap_name, colormap_field) = colormap
                if picked: # every colormap has picked color at max range + 1.0
                    colormap_val = (self.settings["colormaps"][colormap_name][-1][0]+1.0, 0.0, 0.0)
                else:
                    colormap_val = (np.linalg.norm(at.arrays[colormap_field]), 0.0, 0.0)

            for (p, s) in vis_images[i_at]:
                points_lists[at_type].append(p)
                radius_lists[at_type].append(r)
                colormap_vals_lists[at_type].append(colormap_val)
                i_at_lists[at_type].append(i_at)

        return (unique_types, points_lists, radius_lists, colormap_vals_lists, i_at_lists)

    # need to see what can be optimized if settings_only is True
    def update_atom_spheres(self, at):

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
            actor.SetProperty(self.settings["atom_types"][at_type]["prop"])

    # need to see what can be optimized if settings_only is True
    # improve speed by using vtkProgrammableGlyphs like bonds?
    def update_vectors(self, at):
        for actor in self.vector_actors:
            self.renderer.RemoveActor(actor)
        self.vector_actors = []

        if at.info.get("_vtk_vectors", None) is None:
            return

        pos = at.get_positions()
        cell = at.get_cell()
        cell_inv = at.get_reciprocal_cell().T
        rv = at.info.get("_vtk_images", (0.0, 1.0, 0.0, 1.0, 0.0, 1.0))

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

        vis_images = self.visible_images(at)

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

            for (p, s) in vis_images[i_at]:
                if orientation is not None:
                    actors = [vtk.vtkFollower(),vtk.vtkFollower(),vtk.vtkFollower()]
                    for actor in actors:
                        actor.SetCamera(self.renderer.GetActiveCamera())
                        actor.SetPosition(p)
                        # assumes source_orient = \hat{x}

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
    def update_bonds(self, at, color_only):
        if color_only: # colors are handled by fixed vtkProperty, no need to update here
            return

        if not hasattr(at, "bonds"):
            for actor in self.bonds_actors.values():
                actor.VisibilityOff()
                self.renderer.RemoveActor(actor)
            return

        cell = at.get_cell()
        pos = at.get_positions()
        cell_inv = at.get_reciprocal_cell().T
        rv = at.info.get("_vtk_images", (0.0, 1.0, 0.0, 1.0, 0.0, 1.0))

        vis_images = self.visible_images(at)

        vis_images_p = [ [] for i in range(len(at)) ]
        vis_images_s = [ [] for i in range(len(at)) ]
        for i_at in range(len(at)):
            for (p, s) in vis_images[i_at]:
                vis_images_p[i_at].append(p)
                vis_images_s[i_at].append(s)

        points_lists = {}
        axes_lists = {}
        angles_lists = {}
        scales_lists = {}
        i_at_bonds = {}
        glyphs = {}
        for name in self.bond_prop:
            points_lists[name] = []
            axes_lists[name] = []
            angles_lists[name] = []
            scales_lists[name] = []
            i_at_bonds[name] = []
            points_lists[name+"_picked"] = []
            axes_lists[name+"_picked"] = []
            angles_lists[name+"_picked"] = []
            scales_lists[name+"_picked"] = []
            i_at_bonds[name+"_picked"] = []

        for i_at in range(len(at)):
            for (i_bond, b) in enumerate(at.bonds[i_at]):
                j_at = b["j"]
                (dr, dr_norm) = bond_vector(cell, pos, i_at, j_at, b["S"], dist=True)
                name = b["name"]
                rad = self.bond_prop[name].radius

                if b["picked"]:
                    name += "_picked"

                axis = np.cross(dr, [0.0, 1.0, 0.0])
                if np.linalg.norm(axis) < 1.0e-4: # essentially parallel to \hat{y}
                    axis = (1.0,0.0,0.0)
                    angle = 0.0
                else:
                    axis /= np.linalg.norm(axis)
                    dr_hat = dr / dr_norm
                    # angle = -np.arccos(np.dot(dr_hat, [0.0, 1.0, 0.0]))*180.0/np.pi
                    angle = -np.arccos(dr_hat[1])*180.0/np.pi

                for (pi, si) in zip(vis_images_p[i_at], vis_images_s[i_at]):
                    if i_at <= j_at or tuple(si+b["S"]) not in vis_images_s[j_at]:
                        points_lists[name].append(pi+dr/4.0)
                        axes_lists[name].append(axis)
                        angles_lists[name].append(angle)
                        scales_lists[name].append((rad, dr_norm/2.0, rad))
                        i_at_bonds[name].append((i_at, i_bond))

                for (pj, sj) in zip(vis_images_p[j_at], vis_images_s[j_at]):
                    if i_at <= j_at or tuple(sj-b["S"]) not in vis_images_s[i_at]:
                        points_lists[name].append(pj-dr/4.0)
                        axes_lists[name].append(axis)
                        angles_lists[name].append(angle)
                        scales_lists[name].append((rad, dr_norm/2.0, rad))
                        i_at_bonds[name].append((i_at, i_bond))

        # remove unused
        for name in self.bonds_actors:
            if name not in points_lists or len(points_lists[name]) == 0:
                self.bonds_actors[name].VisibilityOff()
                self.renderer.RemoveActor(self.bonds_actors[name])

        for name in points_lists:
            points_list = points_lists[name]
            if len(points_list) == 0:
                continue
            axes_list = axes_lists[name]
            angles_list = angles_lists[name]
            scales_list = scales_lists[name]
            i_at_bond = i_at_bonds[name]
            if name not in self.bonds_points_data:
                self.bonds_points_data[name] = vtk.vtkPoints()
            points = self.bonds_points_data[name]

            if name not in self.bonds_actors:
                self.bonds_actors[name] = vtk.vtkActor()
                actor = self.bonds_actors[name]
                actor._vtk_type = "bonds_glyphs"
                actor.point_to_input_point = self.shapes["sources"]["cylinder_y"][1]
                actor.PickableOn()
            else:
                actor = self.bonds_actors[name]

            # data needed to place, rotate, and scale glyphs
            points.Reset()
            points.SetNumberOfPoints(len(points_list))
            for pt_i in range(len(points_list)):
                points.SetPoint(pt_i, points_list[pt_i])
            axes = numpy_to_vtk(np.array(axes_list))
            angles = numpy_to_vtk(np.array(angles_list))
            scales = numpy_to_vtk(np.array(scales_list))

            glyphs_data = vtk.vtkPolyData()
            glyphs_data.SetPoints(points)
            glyphs_data.GetPointData().AddArray(axes)
            glyphs_data.GetPointData().AddArray(angles)
            glyphs_data.GetPointData().AddArray(scales)

            glyphs = vtk.vtkProgrammableGlyphFilter()
            glyphs.SetSourceConnection(self.shapes["sources"]["cylinder_y"][0].GetOutputPort())
            glyphs.SetInputData(glyphs_data)

            def place_glyph(glyphs=glyphs): # avoid late binding
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

            glyphs.SetGlyphMethod(place_glyph)

            id_glyphs = vtk.vtkIdFilter()
            id_glyphs.SetInputConnection(glyphs.GetOutputPort())
            id_glyphs.SetIdsArrayName("IDs")
            id_glyphs.FieldDataOn()

            glyphs_mapper = vtk.vtkPolyDataMapper()
            glyphs_mapper.SetInputConnection(id_glyphs.GetOutputPort())

            actor.SetMapper(glyphs_mapper)
            actor.i_at_bond = i_at_bond

            if name.endswith("_picked"):
                actor.SetProperty(self.settings["picked"]["prop"])
            else:
                actor.SetProperty(self.bond_prop[name])

            actor.VisibilityOn()
            self.renderer.AddActor(actor)

    def update_legend(self, at):

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
            self.legend_label_actors.append(vtk.vtkActor2D())
            self.legend_label_actors[-1]._vtk_type = "legend_label"

        sphere_scale = self.settings["legend"]["sphere_scale"]

        # create all actors, except for positions, to get actual screen sizes
        sphere_max_size = 0.0
        legend_max_width = 0.0
        legend_max_height = 0.0
        for (l_i, (legend_sphere_actor, legend_label_actor, at_type)) in enumerate(zip(self.legend_sphere_actors, self.legend_label_actors, unique_atom_types)):
            atom_i_of_type = atom_type_list.index(at_type)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(self.shapes["sources"]["sphere"][0].GetOutputPort())
            legend_sphere_actor.SetMapper(mapper)

            prop = get_atom_prop(self.settings, at_type, atom_i_of_type, at.arrays)

            r = sphere_scale*get_atom_radius(self.settings, at_type, list(atom_type_list), at)
            legend_sphere_actor.SetProperty(prop)
            legend_sphere_actor.SetScale(r, r, r)

            legend_sphere_actor.PickableOff()
            legend_sphere_actor.VisibilityOn()

            m = vtk.vtkTextMapper()
            m.SetInput(at_type)
            m.SetTextProperty(self.settings["atom_label"]["prop"])
            legend_label_actor.SetMapper(m)
            legend_label_actor.PickableOff()
            legend_label_actor.VisibilityOn()

        # figure out max sizes
        for actor in self.legend_sphere_actors + self.legend_label_actors:
            self.renderer.AddActor(actor)
        self.renderer.GetRenderWindow().Render()
        max_sphere_size = 0.0
        max_text_size_x = 0.0
        max_text_size_y = 0.0
        for (l_i, (legend_sphere_actor, legend_label_actor, at_type)) in enumerate(zip(self.legend_sphere_actors, self.legend_label_actors, unique_atom_types)):
            b = np.array(legend_sphere_actor.GetBounds())
            max_sphere_size = max(max_sphere_size, np.linalg.norm(b[1::2]-b[0::2])/np.sqrt(3))
            max_text_size_x = max(max_text_size_x, legend_label_actor.GetMapper().GetWidth(self.renderer)*dp_world)
            max_text_size_y = max(max_text_size_y, legend_label_actor.GetMapper().GetHeight(self.renderer)*dp_world)

        # position and space legend using actual sizes
        legend_raw_pos = self.settings["legend"]['position']
        spacing = self.settings["legend"]["spacing"] * max(max_sphere_size*1.5, max_text_size_y*1.5)/dp_world
        legend_offset = np.zeros((2),dtype=int)
        if legend_raw_pos[0] >= 0:
            legend_offset[0] = int(1.3 * max_sphere_size/dp_world)
        else:
            legend_offset[0] = -int( 1.3 * max_sphere_size/dp_world + max_text_size_x/dp_world)
        if legend_raw_pos[1] >= 0:
            legend_offset[1] = int(spacing * (len(unique_atom_types)-0.5))
        else:
            legend_offset[1] = -int(0.5*spacing)
        legend_pos_display = (legend_raw_pos + legend_offset) % display_size

        label_x_offset = int(max_sphere_size/dp_world)
        label_y_offset = -int(max_text_size_y/2.0/dp_world)

        for (l_i, (legend_sphere_actor, legend_label_actor, at_type)) in enumerate(zip(self.legend_sphere_actors, self.legend_label_actors, unique_atom_types)):
            # let's hope atom 0 isn't getting clipped
            self.renderer.SetWorldPoint(list(pos[0]) + [1.0])
            self.renderer.WorldToDisplay()
            arb_point = self.renderer.GetDisplayPoint()

            # got from legend pos in display coords to world pos (for sphere)
            sphere_pos_display = [legend_pos_display[0], legend_pos_display[1]-spacing*l_i, arb_point[2]]
            self.renderer.SetDisplayPoint(sphere_pos_display)
            self.renderer.DisplayToWorld()
            sphere_pos_world = self.renderer.GetWorldPoint()

            legend_sphere_actor.SetPosition(np.array(sphere_pos_world[0:3])/sphere_pos_world[3])

            label_pos_display = [ sphere_pos_display[0]+label_x_offset, sphere_pos_display[1]+label_y_offset]
            legend_label_actor.SetPosition(label_pos_display)


    def string_dollar_sub(self, string, at, i_at=None):
        dict_frame = at.info
        dict_frame_special_case = { "config_n" : self.cur_frame }
        dict_atom = at.arrays
        try:
            magmoms = at.get_magnetic_moments()
        except:
            magmoms = None
        try:
            init_magmoms = at.get_initial_magnetic_moments()
        except:
            init_magmoms = None
        dict_atom_special_case = { "ID" : list(range(len(at))), 
                                   "Z" : at.get_atomic_numbers(),
                                   "species" : at.get_chemical_symbols(),
                                   "magmom" : magmoms,
                                   "initial_magmoms" : init_magmoms
                                 }
        # substitute $${} from at.arrays
        for substr in re.findall('\$\${([^}]*)}', string):
            if i_at is None:
                raise ValueError("found $${{}} substitution for '{}' but i_at is None. Perhaps incorrectly set a per-atom property ($${...}) in a per-frame context?".format(substr))
            if dict_atom_special_case is not None and substr in dict_atom_special_case:
                string = re.sub('\$\${'+substr+'}',str(dict_atom_special_case[substr][i_at]),string,count=1)
            elif dict_atom is not None and substr in dict_atom:
                string = re.sub('\$\${'+substr+'}',str(dict_atom[substr][i_at]),string,count=1)
            else:
                raise ValueError("Tried to do $${{}} substitution on undefined per-atom field '{}'".format(substr))
        # substitute ${} from at.info
        for substr in re.findall('(?:^|[^$])\${([^}]*)}', string):
            if dict_frame_special_case is not None and substr in dict_frame_special_case:
                string = re.sub('\${'+substr+'}',str(dict_frame_special_case[substr]),string,count=1)
            elif dict_frame is not None and substr in dict_frame:
                string = re.sub('\${'+substr+'}',str(dict_frame[substr]),string,count=1)
            else:
                raise ValueError("Tried to do ${{}} substitution on undefined per-frame field '{}'".format(substr))

        # eval $()
        while '$(' in string:
            p0 = string.index('$(')
            depth = 0
            for char_i in range(p0+2, len(string)):
                if string[char_i] == '(':
                    depth += 1
                elif string[char_i] == ')':
                    if depth > 0:
                        depth -= 1
                    else: # depth == 0
                        break
            if char_i == len(string)-1 and depth > 1:
                raise RuntimeError("Failed to find closing of '$(' in '{}'".format(string))
            substring = string[p0+2:char_i]
            string = string.replace('$('+substring+')',str(eval(substring)))

        return string

    def update_frame_label(self, at):
        if self.frame_label_actor is None:
            self.frame_label_actor = vtk.vtkTextActor()
            self.frame_label_actor._vtk_type = "frame_label"
            self.frame_label_actor.PickableOff()

        if not self.settings["frame_label"]["show"]:
            self.frame_label_actor.VisibilityOff()
            return

        self.frame_label_actor.VisibilityOn()

        self.frame_label_actor.SetTextProperty(self.settings["frame_label"]["prop"])
        self.frame_label_actor.SetDisplayPosition(20,20)

        if at.info.get("_vtk_frame_label_string", None) is not None:
            frame_label_string = at.info["_vtk_frame_label_string"]
        else:
            frame_label_string = self.settings["frame_label"]["string"]

        if frame_label_string == "_NONE_":
            frame_label_str = ""
        else:
            frame_label_str = self.string_dollar_sub(frame_label_string, at)

        self.frame_label_actor.SetInput(frame_label_str)
        self.renderer.AddActor(self.frame_label_actor)

    def vector_lut(self, colors):
        # look up table for bond colors
        vector_lut = vtk.vtkColorTransferFunction()
        vector_lut.SetRange(0, len(colors)-1)
        for (i, c) in enumerate(colors):
            vector_lut.AddRGBPoint(float(i), c[0], c[1], c[2])
        return vector_lut

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
            lut.AddRGBPoint(data[-1][0]+1.0, self.settings["picked"]["color"][0], self.settings["picked"]["color"][1], self.settings["picked"]["color"][2])
            self.luts[name] = lut

        self.atom_type_luts = {}
        # create trivial color maps for fixed colors, or link to true maps
        for name in self.settings["atom_types"].get_all():
            color = self.settings["atom_types"][name]["color"]
            if color is not None: # just a color
                lut = vtk.vtkColorTransferFunction()
                lut.SetRange(0,0)
                lut.AddRGBPoint(0.0, color[0], color[1], color[2])
                lut.AddRGBPoint(1.0, self.settings["picked"]["color"][0], self.settings["picked"]["color"][1], self.settings["picked"]["color"][2])
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
    def update_atom_labels(self, at):
        for actor in self.atom_label_actors:
            self.renderer.RemoveActor(actor)
        self.atom_label_actors = []

        if "_vtk_atom_label_show" in at.info and at.info["_vtk_atom_label_show"] is not None:
            show_labels = at.info["_vtk_atom_label_show"]
        else:
            show_labels = self.settings["atom_label"]["show"]
        if show_labels is False:
            return

        atom_type_list = get_atom_type_list(self.settings, at)
        vis_images = self.visible_images(at)

        dp_world = self.label_offset_world(at.positions[0])

        IDs = list(range(len(at)))
        Zs = at.get_atomic_numbers()
        species = at.get_chemical_symbols()
        for i_at in IDs:
            label_raw_string = None
            if "_vtk_label" in at.arrays: # try per-atom value first
                label_raw_string = at.arrays["_vtk_label"][i_at]
                if re.search('^\s*$', label_raw_string) or re.search('^"?\s*"?$', label_raw_string) or re.search("^'?\s*'?$", label_raw_string):
                    label_raw_string = None
            if label_raw_string is None: # use string from per-config or global setting
                if "_vtk_atom_label_string" in at.info:
                    label_raw_string = at.info["_vtk_atom_label_string"]
                else:
                    label_raw_string = self.settings["atom_label"]["string"]
            if label_raw_string == "_NONE_":
                label_str = ""
            else:
                label_str = self.string_dollar_sub(label_raw_string, at, i_at)
            
            r = get_atom_radius(self.settings, atom_type_list[i_at], i_at, at)
            if dp_world > 0:
                dp_disp = 3+0.7*r/dp_world
            else:
                dp_disp = 0

            for (p, s) in vis_images[i_at]:
                label_actor = vtk.vtkBillboardTextActor3D()
                label_actor._vtk_type="atom_label"
                label_actor.PickableOff()
                label_actor.SetVisibility(True)

                label_actor.SetInput(label_str)
                label_actor.SetPosition(p)
                label_actor.SetDisplayOffset(int(dp_disp), int(dp_disp))
                label_actor.SetTextProperty(self.settings["atom_label"]["prop"])

                self.renderer.AddActor(label_actor)
                self.atom_label_actors.append(label_actor)

    def measure(self, n=None, frame_i=None):
        if frame_i is None:
            at = self.at_list[self.cur_frame]
        else:
            at = self.at_list[frame_i]

        if n is None:
            if "_vtk_picked" in at.arrays:
                at_indices = np.where(at.arrays["_vtk_picked"])[0]
            else:
                at_indices = []
        else:
            at_indices = n

        p = at.get_positions()
        c = at.get_cell()

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
                        (dr, dist) = bond_vector(c, p, i_at, b["j"], b["S"], dist=True)
                        print("bond {} {} vec {} {} {} ({})".format(i_at, b["j"], dr[0], dr[1], dr[2], dist))
            for i_at in range(len(at)):
                for b1 in at.bonds[i_at]:
                    j_at = b1["j"]
                    if b1["picked"]:
                        for b2 in at.bonds[i_at]:
                            k_at = b2["j"]
                            if b2["picked"] and k_at <= j_at and (k_at != j_at or np.any(b1["S"] != -b2["S"])):
                                (v1, v1_dist) = bond_vector(c, p, i_at, b1["j"], b1["S"], dist=True)
                                (v2, v2_dist) = bond_vector(c, p, i_at, b2["j"], b2["S"], dist=True)
                                ang = 180.0/np.pi * np.arccos(np.dot(v1,v2)/(v1_dist*v2_dist))
                                print("bond-angle {} {} {} angle {}".format(i_at, j_at, k_at, ang))

    def supercell(self, n_dup, wrap=True, frames=None):
        for frame_i in self.frame_list(frames):
            at = self.at_list[frame_i]

            # store original info if this is first duplicate
            if "dup_orig_n" not in at.info:
                at.info["dup_orig_n"] = len(at)
            if "dup_orig_cell" not in at.info:
                at.info["dup_orig_cell"] = at.get_cell()
            if "dup_orig_index" not in at.arrays:
                at.new_array("dup_orig_index",np.array(range(len(at))))

            if len(n_dup) == 3:
                if np.any(np.logical_and(np.logical_not(at.get_pbc()),n_dup != 1)):
                    raise RuntimeError("Trying to duplicate along direction that is not periodic")
                at *= n_dup
            else:
                P = np.reshape(n_dup,(3,3))
                if np.dot(P[0],np.cross(P[1],P[2])) < 0:
                    raise RuntimeError("Supercell vectors need have positive scalar triple product, reorder 2 or flip 1")
                new_at = make_supercell(at, P)

                # eliminate atoms in new_at that are at same positions as old at
                cinv = new_at.get_reciprocal_cell().T
                at_in_new_lat = np.matmul(at.get_positions(), cinv)
                new_lat_pos = new_at.get_scaled_positions()
                keep_inds = list(range(len(new_at)))
                for i_at in range(len(at)):
                    d = new_lat_pos - at_in_new_lat[i_at]
                    if np.min(np.linalg.norm(d-np.round(d),axis=1)) < 1.0e-6:
                        keep_inds.remove(i_at)

                at.set_cell(new_at.get_cell(), False)
                at += new_at[keep_inds]

            if wrap:
                at.wrap()

            # maybe some day figure out how to duplicate bonds
            # for now just get rid of all of them
            if hasattr(at, "bonds"):
                at.bonds.reinit()

        self.update(frames)

    def bond(self, at, name, at_type1, at_type2, criterion, no_pbc=False):
        if name is None:
            name = "default"

        if not hasattr(at, "bonds"):
            at.bonds = DavTKBonds(at, self.settings)

        if criterion == "picked":
            if "_vtk_picked" in at.arrays:
                indices = np.where(self.cur_at().arrays["_vtk_picked"])[0]
            else:
                indices = []
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
            at.bonds.cutoff(name, criterion[1], at_type1, at_type2, across_pbc=not no_pbc)
        else:
            raise ValueError("Unknown bonding criterion type '{}'".format(criterion[0]))

    def snapshot(self, filename=None, mag=1):
        renderLarge = vtk.vtkRenderLargeImage()
        renderLarge.SetInput(self.renderer)
        renderLarge.SetMagnification(mag)

        if filename is not None:
            if not filename.endswith(".png"):
                warn("creating PNG file with suffix that is not .png") 
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

    def delete_volume_rep(self, at, name):
        for (actor, _) in at.volume_reps[name][2]:
            self.renderer.RemoveActor(actor)
        del at.volume_reps[name]

    def make_isosurface(self, img, transform, params, prop_name):

        (threshold,) = params

        isosurface = vtk.vtkMarchingCubes()
        isosurface.SetInputData( img.GetOutput() )
        isosurface.ComputeNormalsOn()
        isosurface.SetValue( 0, threshold )
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
        actor.SetProperty(self.volume_rep_prop[prop_name])
        actor.GetProperty().BackfaceCullingOff()
        actor.PickableOff()

        return actor

    def add_volume_rep(self, name, data, style, params, cmd_string):
        # at.volume_reps[name] is 3-element tuple: image, transform (only one for each data set), 
        #     and list of (actor, command_string) tuples (could be multiple views of same data).
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
            at.volume_reps[name][2].append((self.make_isosurface(img, transform, params, name), cmd_string))
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

            # save vectors
            if "_vtk_vectors" in at.info:
                at.info["__vtk_vectors_INFO"] = True
                for (k,v) in at.info["_vtk_vectors"].items():
                    at.info["_vtk_vectors_"+k] = v

            if "_vtk_commands" not in at.info:
                at.info["_vtk_commands"] = ""

            # save volume commands
            if hasattr(at,"volume_reps") and len(at.volume_reps) > 0:
                for volume_rep in at.volume_reps.values():
                    for (_, cmd) in volume_rep[2]:
                        at.info["_vtk_commands"] += ' '.join(cmd) + " ; "

            # save bond, polyhedra, and volume_rep properties
            for (cmd, prop_dict) in [("bond", self.bond_prop), ("polyhedra", self.polyhedra_prop), ("volume", self.volume_rep_prop)]:
                for (name, prop) in prop_dict.items():
                    col = prop.GetColor()
                    at.info["_vtk_commands"] += "{} -name {} -color {} {} {} -opacity {} -specular {} -specular_radius {} -ambient {}".format(
                        cmd, name, col[0], col[1], col[2], prop.GetOpacity(), prop.GetSpecular(), 1.0/prop.GetSpecularPower(), prop.GetAmbient())
                    try:
                        at.info["_vtk_commands"] += " -r {}".format(prop.radius)
                    except AttributeError:
                        pass
                    at.info["_vtk_commands"] += " ; "

        for at in ats:
            if at.info["_vtk_commands"].endswith("; "):
                at.info["_vtk_commands"] = at.info["_vtk_commands"][:-2]
            if re.search('^[\s;]*$', at.info["_vtk_commands"]):
                del at.info["_vtk_commands"]

    def prep_after_atoms_read(self, ats=None):
        if ats is None:
            ats = self.at_list

        for at in ats:
            # reconstruct bonds from arrays
            if "_vtk_bonds" in at.arrays:
                at.bonds = DavTKBonds(at, self.settings)
                at.bonds.read_from_atoms_arrays()
            # reconstruct vectors from info
            if at.info.get("__vtk_vectors_INFO",False):
                at.info["_vtk_vectors"] = {}
                for (k,v) in at.info.items():
                    if k.startswith("_vtk_vectors_"):
                        at.info["_vtk_vectors"][k.replace("_vtk_vectors_","")] = v
                del at.info["__vtk_vectors_INFO"]

    def arb_polyhedra(self, at, name, atom_lists):
        point_strings = [ "" ] * len(at)
        face_strings = [ "" ] * len(at)
        for l in atom_lists:
            for i in l:
                point_strings[l[0]] += "_".join(["{}".format(x) for x in at.get_distance(l[0], i, mic=True, vector=True)]) + "_"
            face_strings[l[0]] += "_".join(["{}".format(i) for i in l]) + "__"
        for i in range(len(at)):
            point_strings[i] = re.sub(r'_$', '', point_strings[i])
            face_strings[i] = re.sub(r'__$', '', face_strings[i])

        strings = []
        for i in range(len(at)):
            if len(point_strings[i]) > 0 and len(face_strings[i]) > 0:
                strings.append(point_strings[i]+"___"+face_strings[i])
            else:
                strings.append("_NONE_")

        try:
            del at.arrays["_vtk_polyhedra_"+name]
        except KeyError:
            pass
        at.new_array("_vtk_polyhedra_"+name, np.array(strings))

    def coordination_polyhedra(self, at, name, center_at_type, neighb_at_type, cutoff=None, bond_name=None):
        if cutoff is None: # no cutoff, use existing bonds
            if not hasattr(at, "bonds"):
                warn("coordination polyhedra without cutoff require bonds already exist") 
                return
            if bond_name is None:
                bond_name = "default"
        else: # cutoff is specified
            if bond_name is not None:
                raise ValueError("polyhedra got both cutoff {} and bond_name '{}'".format(cutoff,bond_name))

        pos = at.get_positions()
        cell = at.get_cell()
        points = [ [] for i in range(len(at)) ]
        strings = [ "" ] * len(at)
        atom_type_list = get_atom_type_list(self.settings, at)

        if cutoff is None:
            for i_at in range(len(at)):
                if atom_type_list[i_at] != center_at_type:
                    continue
                for b in at.bonds[i_at]:
                    j_at = b["j"]
                    if neighb_at_type is not None and atom_type_list[j_at] != neighb_at_type:
                        continue
                    if bond_name is not None and b["name"] != bond_name:
                        continue

                    D = bond_vector(cell, pos, i_at, j_at, b["S"])
                    points[i_at].append(D)
        else:
            nn_list = ase.neighborlist.neighbor_list('ijdD', at, cutoff, self_interaction=True)
            for (i,j,d,D) in zip(nn_list[0], nn_list[1], nn_list[2], nn_list[3]):
                if atom_type_list[i] != center_at_type or d == 0.0 or (neighb_at_type is not None and atom_type_list[j] != neighb_at_type):
                    continue
                points[i].append(D)

        for i_at in range(len(at)):
            if len(points[i_at]) > 0:
                hull = scipy.spatial.ConvexHull(np.array(points[i_at]))

                for point in points[i_at]:
                    strings[i_at] += "_".join([str(v) for v in point]) + "_"
                strings[i_at] = re.sub(r'_$', '', strings[i_at])

                strings[i_at] += "___"

                for simplex in hull.simplices:
                    strings[i_at] += "_".join([str(i) for i in simplex]) + "__"
                strings[i_at] = re.sub(r'__$', '', strings[i_at])

            else:
                strings[i_at] = "_NONE_"

        # add strings
        try:
            del at.arrays["_vtk_polyhedra_"+name]
        except KeyError:
            pass
        at.new_array("_vtk_polyhedra_"+name, np.array(strings))

    def update_polyhedra(self, at, settings_only = False):
        if settings_only: # colors handled by fixed vtkProperty
            return

        # clean up actors, will be re-added later as needed
        for actor in self.polyhedra_actors:
            self.renderer.RemoveActor(actor)

        if len(self.polyhedra_prop) == 0:
            return

        self.polyhedra_actors = []

        vis_images = self.visible_images(at)
        pos = at.positions

        for name in self.polyhedra_prop:
            if "_vtk_polyhedra_"+name not in at.arrays:
                continue

            for (i_at, string) in enumerate(at.arrays["_vtk_polyhedra_"+name]):
                if string == "_NONE_":
                    continue
                m = re.search('(.*)___(.*)', string)
                (points_str, polygons_str) = (m.group(1), m.group(2))
                points_data = np.reshape([float(p) for p in points_str.split("_")],(-1,3))

                for (p, s) in vis_images[i_at]:
                    points = vtk.vtkPoints()
                    for D in points_data:
                        points.InsertNextPoint(p+D)
                    polygons = vtk.vtkCellArray()
                    for m in polygons_str.split("__"):
                        polygon_data = [int(i) for i in m.split("_")]
                        polygon = vtk.vtkPolygon()
                        polygon.GetPointIds().SetNumberOfIds(len(polygon_data))
                        for (i_vertex, i) in enumerate(polygon_data):
                            polygon.GetPointIds().SetId(i_vertex, i)
                        polygons.InsertNextCell(polygon)

                    polydata = vtk.vtkPolyData()
                    polydata.SetPoints(points)
                    polydata.SetPolys(polygons)

                    mapper = vtk.vtkPolyDataMapper()
                    mapper.SetInputData(polydata)

                    actor = vtk.vtkActor()
                    actor._vtk_type = "polyhedron"
                    actor.SetMapper(mapper)
                    actor.SetProperty(self.polyhedra_prop[name])
                    actor.PickableOff()

                    self.renderer.AddActor(actor)
                    self.polyhedra_actors.append(actor)

    def set_view(self, along_up, lattice, scale):
        cam = self.renderer.GetActiveCamera()

        if along_up is not None:
            along = along_up[0:3]
            up = along_up[3:6]

            if lattice:
                cell = self.cur_at().get_cell()
                along = np.dot(along, cell)
                up = np.dot(up, cell)
            along /= np.linalg.norm(along)
            up /= np.linalg.norm(up)

            focal_point = np.array(cam.GetFocalPoint())
            pos = np.array(cam.GetPosition())
            dist = np.linalg.norm(pos-focal_point)
            new_pos = focal_point - dist * along

            cam.SetPosition(new_pos)
            cam.SetViewUp(up)

        if scale is not None:
            cam.SetParallelScale(cam.GetParallelScale()/scale)
