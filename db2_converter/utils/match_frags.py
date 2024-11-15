from rdkit import Chem
from rdkit.Chem.rdMolAlign import AlignMolConformers

from db2_converter.utils.utils import next_mol2_lines
sb_smarts = "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"


def getHeavyAtomNeighbors(atom1):
    return [n for n in atom1.GetNeighbors() if n.GetSymbol() != "H"]


def find_central_atom(mol, distmat):
    sums = []
    if mol.GetNumHeavyAtoms() <= 2:
        for i in range(mol.GetNumAtoms()):
            # only consider non-terminal atoms
            if mol.GetAtomWithIdx(i).GetSymbol() != "H":
                aid1=i
    else:
        for i in range(mol.GetNumAtoms()):
            # only consider non-terminal atoms
            if len(getHeavyAtomNeighbors(mol.GetAtomWithIdx(i))) < 2:
                continue
            tmp = [d for d in distmat[i]]
            tmp.pop(i)
            sums.append((sum(tmp), i))
        sums.sort()
        aid1 = sums[0][1]
    return aid1  # most central atom comes first


def find_central_rigid(mol):
    distmat = Chem.GetDistanceMatrix(mol)
    aid1 = find_central_atom(mol, distmat)
    aid_others = [
        atom.GetIdx() for atom in getHeavyAtomNeighbors(mol.GetAtomWithIdx(aid1))
    ]
    return [aid1, *aid_others]


def dist_2(c1, c2):
    """Calculate the distance**2 between two points"""
    return (c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2


def get_same_index(ref, i1, i2, index, eps, eps_modify):
    """A function that checks whether two conformers share the same coordinates for a part of the whole molecule.
    Here, two cutoffs are applied:
    EPS: The cutoff to define two coordinate are different.
    EPS_modify: two coordinate distance > EPS_modify are also regarded as same. This script will assign the exact same coordinates to the second conformer to the first conformer. Thus, two conformer would share the same coordinates.
    Other arguments:
    ref:
    """
    same_index = list()
    not_same = list()
    eps_2 = eps**2
    eps_modify_2 = eps_modify**2
    conf1 = ref.GetConformer(id=i1)
    conf2 = ref.GetConformer(id=i2)
    for i in index:
        tmp = conf1.GetAtomPosition(i)
        c1 = tmp.x, tmp.y, tmp.z
        tmp = conf2.GetAtomPosition(i)
        c2 = tmp.x, tmp.y, tmp.z
        if dist_2(c1, c2) > eps_2:
            not_same.append(i)
        elif eps_modify_2 < dist_2(c1, c2) <= eps_2:
            same_index.append(
                i
            )  # see them as same coordinate. But new conformation should be assign IDENTICAL coordinate.
            conf2.SetAtomPosition(i, conf1.GetAtomPosition(i))  # Assign
        else:
            same_index.append(i)

    if len(same_index) <= 1:
        return False, same_index
    else:  # 2 atoms matched is enough for is_same
        return True, same_index


def f_AlignMolConformers(ref, atomIds, eps, eps_modify, maxIters=500, reflect=False):
    """Align conformers using rdkit.Chem.rdMolAlign.AlignMolConformers"""
    index = atomIds
    rmsd = list()
    if (
        len(ref.GetConformers()) == 1
    ):  # if Only one conformation, then no alignation is needed.
        return [[0]]
    else:
        # RMSlist if provided, fills in the RMS values between the reference
        # conformation and the other aligned conformations
        # RMSlist will be written in rmsd
        AlignMolConformers(
            ref, atomIds=atomIds, maxIters=maxIters, reflect=reflect, RMSlist=rmsd
        )  # After this func, rmsd should be a list with length=conf_num-1
        all_families = list()
        all_indexes = list()
        # Go through all conformations
        for i in range(len(rmsd) + 1):  # len(rmsd)+1 is exactly the conf_num
            flag_match = False
            # Go match all exist families ( But match most one )
            for current_family, current_index, tmpi in zip(
                all_families, all_indexes, range(len(all_families))
            ):  # only 3 items are all non-null, the for-loop can go below
                fam_ref = current_family[0]
                fam_index = current_index
                same_flag, same_index = get_same_index(ref, fam_ref, i, eps=eps, eps_modify=eps_modify, index=fam_index)
                if same_flag and len(same_index) == len(
                    index
                ):  # But here, only all atoms in index are matched are considered as True
                    flag_match = True
                    all_families[tmpi].append(
                        i
                    )  # if flag_matched, all_families will have element list with len>1
                    break
            if not flag_match:  # Do not match, Create new family
                all_families.append([i])
                all_indexes.append(atomIds[:])  # each fragindex (a group of indexes)
        return all_families


def create_editable_mol(mol, bonds):
    em = Chem.EditableMol(mol)
    nAts = mol.GetNumAtoms()
    for a, b in bonds:
        em.RemoveBond(a, b)
        em.AddAtom(Chem.Atom(0))
        em.AddBond(a, nAts, Chem.BondType.SINGLE)
        em.AddAtom(Chem.Atom(0))
        em.AddBond(b, nAts + 1, Chem.BondType.SINGLE)
        nAts += 2
    fmol = em.GetMol()
    Chem.SanitizeMol(fmol)
    return fmol


def embed_blocks_molconformer(blocks, removeHs=False):
    ref_mol = Chem.MolFromMol2Block("".join(blocks[0]), removeHs=removeHs)
    if len(blocks) == 1:
        return ref_mol
    all_mols = [Chem.MolFromMol2Block("".join(x), removeHs=False) for x in blocks[1:]]
    # add all molecules to ref molecule's conformations
    for mol in all_mols:
        mol_conf = mol.GetConformer()
        ref_mol.AddConformer(mol_conf, assignId=True)
    return ref_mol

def mol_to_frags(mol, cutsmarts=sb_smarts):
    patt = Chem.MolFromSmarts(cutsmarts)  # acyclic single bond
    bonds = mol.GetSubstructMatches(patt)
    # create an editable molecule, break the bonds, and add dummies:
    fmol = create_editable_mol(mol, bonds)
    frags = [x for x in Chem.GetMolFrags(fmol, asMols=True)]
    return frags

def mol_to_ring_frags(mol, cutsmarts=sb_smarts, minfragsize=3):
    frags = mol_to_frags(mol, cutsmarts=cutsmarts)
    suitable_ring_frags = []
    for frag in frags:
        if not (frag.GetRingInfo().AtomRings() and (len(frag.GetAtoms()) >= minfragsize)): # frag must be ring >= 3 heavy atoms
            continue
        else:
            suitable_ring_frags.append(frag)
    return suitable_ring_frags

def uniquify_hyd_rotamers_mol2file(inmol2,outmol2):
    all_blocks = [x for x in next_mol2_lines(inmol2)]
    mol = Chem.MolFromMol2Block("".join(all_blocks[0]), removeHs=False)
    mol_withconfs = embed_blocks_molconformer(all_blocks, removeHs=False)
    heavy_atom_index = [ atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1 ]
    EPS = 0.05
    EPS_modify = 1e-4
    all_index = f_AlignMolConformers(
        mol_withconfs,
        atomIds=heavy_atom_index,
        maxIters=500,
        reflect=False,
        eps=EPS,
        eps_modify=EPS_modify,
    )
    unique_block_indexes = [ index[0] for index in all_index ]
    unique_blocks = [ all_blocks[i] for i in unique_block_indexes ]
    with open(outmol2, "w") as f:
        for block in unique_blocks:
            f.write("".join(block))
    return len(all_blocks) % len(unique_blocks) == 0 # If not zero, means something unexpected happened in this process
