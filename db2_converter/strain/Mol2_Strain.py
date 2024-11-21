from db2_converter.mol2db2_py3_strain import mol2
from db2_converter.strain.Torsion_Strain import calc_strain
from db2_converter.strain.TL_Functions import Mol2MolSupplier

def mol2_strain(inputmol2file):
    db2in_strain = Mol2MolSupplier(inputmol2file) # names and mols
    db2in_standard = mol2.Mol2(inputmol2file)
    tE, pE = calc_strain(*db2in_strain)
    db2in_standard.addStrainInfo(tE, pE)
    return db2in_standard
