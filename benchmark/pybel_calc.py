from openbabel import pybel
import sys
import time

start = time.process_time()

mol = next(pybel.readfile("pdb", sys.argv[1]))
mol.OBMol.PerceiveBondOrders()
ff = pybel._forcefields["gaff"]
if not ff.Setup(mol.OBMol):
    raise RuntimeError("GAFF setup failed for this molecule")

energy = ff.Energy()   # in kcal/mol
print("Single-point GAFF energy:", energy * 4.184, "kJ/mol")
end = time.process_time()

print("CPU time:", end - start, "seconds")
