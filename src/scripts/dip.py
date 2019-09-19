import pybel
import openbabel
import math

poly_smiles = 'C1=[S]C(=C(O1)[N]C=O)C=Cc1ccc2c(c1)oc1c2ccc(c1)C1=[S]C(=C(O1)[N]C=O)C=CC1=[S]C(=C(O1)[N]C=O)C=CC1=[S]C(=C(O1)[N]C=O)C=C'

ob = pybel.ob

def make3D(mol):
    '''
    Makes the mol object from SMILES 3D

    Parameters
    ---------
    mol: object
        pybel molecule object
    '''
    pybel._builder.Build(mol.OBMol)
    mol.addh()


mol = pybel.readstring("smi", poly_smiles)
make3D(mol)
cm = ob.OBChargeModel_FindType('mmff94')
cm.ComputeCharges(mol.OBMol)
diptensor = cm.GetDipoleMoment(mol.OBMol)
dipole = math.sqrt(diptensor.GetX()**2 +
       diptensor.GetY()**2 + diptensor.GetZ()**2)

print(dipole)

def globalopt(mol, debug=False, fast=False):
    # Function which does the geometry optimization
    pybel._builder.Build(mol.OBMol)
    mol.addh()

    ff = pybel._forcefields["mmff94"]
    # ff = pybel._forcefields["uff"]
    success = ff.Setup(mol.OBMol)
    if not success:
        ff = pybel._forcefields["uff"]
        success = ff.Setup(mol.OBMol)
        if not success:
            sys.exit("Cannot set up forcefield")

    ff.SteepestDescent(1000, 1.0e-4)
    ff.WeightedRotorSearch(250, 25)
    ff.WeightedRotorSearch(250, 25)
    ff.ConjugateGradients(500, 1.0e-6)
    ff.GetCoordinates(mol.OBMol)


mol = pybel.readstring("smi", poly_smiles)
globalopt(mol)
cm = ob.OBChargeModel_FindType('mmff94')
cm.ComputeCharges(mol.OBMol)
diptensor = cm.GetDipoleMoment(mol.OBMol)
dipole = math.sqrt(diptensor.GetX()**2 +
       diptensor.GetY()**2 + diptensor.GetZ()**2)

print(dipole)
