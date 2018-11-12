import ase, ase.io, sys, numpy as np

ats = ase.io.read("t.xyz",":")
ats[0].arrays["test"] = np.array([False] * len(ats[0]))

ase.io.write(sys.stdout,ats[0],format="extxyz")

del ats[0][2]

ase.io.write(sys.stdout,ats[0],format="extxyz")
