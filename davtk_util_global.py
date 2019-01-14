import numpy as np

def undo_offset(atoms):
    if "_vtk_offset" in atoms.info:
        atoms.set_positions(atoms.get_positions() - atoms.info["_vtk_offset"])
        atoms.info["_vtk_offset"][:] = 0.0

def center_around(atoms, ctr):
    p = atoms.get_positions()
    if isinstance(ctr, int):
        ctr = p[ctr]
    offset = - np.array(ctr) + np.sum(atoms.get_cell(),axis=0)/2.0
    atoms.set_positions(p + offset)
    if "_vtk_offset" not in atoms.info:
        atoms.info["_vtk_offset"] = np.array([0.0]*3)
    atoms.info["_vtk_offset"] += offset
    atoms.wrap()
