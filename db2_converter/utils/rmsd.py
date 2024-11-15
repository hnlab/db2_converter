import math
from rdkit import Chem
from rdkit.Chem.AllChem import AlignMol

from db2_converter.utils.utils import next_mol2_lines


def MolsFromMol2File(mol2file, removeHs=True):
    try:
        mol2_blocks = [x for x in list(next_mol2_lines(mol2file)) if x]
        mol2_mols = [
            Chem.MolFromMol2Block("".join(x), removeHs=removeHs) for x in mol2_blocks
        ]
        return mol2_blocks, mol2_mols
    except:
        return [], []


def GetMaps(refmol, probemol):
    matches = refmol.GetSubstructMatches(
        probemol, uniquify=False
    )  # iterate atom orders
    if not matches:
        raise ValueError(
            f'mol {refmol.GetProp("_Name")} does not match mol {probemol.GetProp("_Name")}'
        )
    if len(matches) > 1e6:
        print(
            f'{len(matches)} matches detected for molecule {refmol.GetProp("_Name")}, this may lead to a performance slowdown.'
        )

    maps = [list(enumerate(match)) for match in matches]
    return maps


def GetBestRMSD(
    probe,
    ref,
    refConfId=-1,
    probeConfId=-1,
    maps=None,
    align=True,
    removeHs=True,
):
    """Returns the optimal RMS for aligning two molecules, taking
    symmetry into account. As a side-effect, the probe molecule is
    left in the aligned state.

    Arguments:
      - ref: the reference molecule
      - probe: the molecule to be aligned to the reference
      - refConfId: (optional) reference conformation to use
      - probeConfId: (optional) probe conformation to use
      - maps: (optional) a list of lists of (probeAtomId,refAtomId)
        tuples with the atom-atom mappings of the two molecules.
        If not provided, these will be generated using a substructure
        search.

    Returns: #modified by qcxia
      - rmsd: calculated rmsd
      - probe: transformed probe mol

    Note:
    This function will attempt to align all permutations of matching atom
    orders in both molecules, for some molecules it will lead to 'combinatorial
    explosion' especially if hydrogens are present.
    Use 'rdkit.Chem.AllChem.AlignMol' to align molecules without changing the
    atom order.

    """
    # Remove Hs of input molecules
    if removeHs:
        probe = Chem.rdmolops.RemoveHs(probe)
        ref = Chem.rdmolops.RemoveHs(ref)

    assert probe

    # When mapping the coordinate of probe will changed!!!
    ref.pos = orginXYZ(ref)
    probe.pos = orginXYZ(probe)
    if not maps:
        maps = GetMaps(ref, probe)

    bestRMSD = 10000.0
    for amap in maps:
        if align:
            rmsd = AlignMol(probe, ref, probeConfId, refConfId, atomMap=amap)
        else:
            rmsd = RMSD_NotAlign(probe, ref, amap)
        bestRMSD = min(bestRMSD, rmsd)
    return bestRMSD


# Map is probe -> ref
# [(1:3),(2:5),...,(10,1)]
def RMSD_NotAlign(probe, ref, amap):
    rmsd = 0.0
    # print(amap)
    atomNum = ref.GetNumAtoms() + 0.0
    for (pi, ri) in amap:
        posp = probe.pos[pi]
        posf = ref.pos[ri]
        rmsd += dist_2(posp, posf)
    rmsd = math.sqrt(rmsd / atomNum)
    return rmsd


def dist_2(atoma_xyz, atomb_xyz):
    dis2 = 0.0
    for i, j in zip(atoma_xyz, atomb_xyz):
        dis2 += (i - j) ** 2
    return dis2


def orginXYZ(mol):
    mol_pos = {}
    for i in range(0, mol.GetNumAtoms()):
        pos = mol.GetConformer().GetAtomPosition(i)
        mol_pos[i] = pos
    return mol_pos


def RMSDfilter(zinc, samplopts, RMSthres):
    for k, samplopt in enumerate(samplopts):
        try:
            ref_mol2_blocks, ref_mol2_mols = MolsFromMol2File(
                f"conformer.{zinc}.{samplopts[k]}.fixed.mol2", removeHs=False
            )
            if ref_mol2_blocks:
                break
        except:
            continue
    out_mol2_mols, out_mol2_blocks = [], []
    for i, moli in enumerate(ref_mol2_mols):
        try:
            min_rmsd = 1e6
            for j in range(i+1, len(ref_mol2_mols)):
                min_rmsd = min(
                    GetBestRMSD(moli, ref_mol2_mols[j], align=True, removeHs=False),
                    min_rmsd,
                )
            if min_rmsd > RMSthres:
                out_mol2_mols.append(moli)
                out_mol2_blocks.append(ref_mol2_blocks[i])
        except:
            continue
    out_mol2_blocks.append(ref_mol2_blocks[-1])
    if len(samplopts) == 1:
        return out_mol2_blocks
    else:
        ref_mol2_mols = out_mol2_mols
        for samplopt in samplopts[k+1:]:
            try:
                probe_mol2_blocks, probe_mol2_mols = MolsFromMol2File(
                    f"conformer.{zinc}.{samplopt}.fixed.mol2", removeHs=False
                )
                if not probe_mol2_mols:
                    continue
                for i, probe_mol in enumerate(probe_mol2_mols):
                    min_rmsd = 1e6
                    for ref_mol in ref_mol2_mols:
                        min_rmsd = min(
                            GetBestRMSD(probe_mol, ref_mol, align=True, removeHs=False),
                            min_rmsd,
                        )
                    if min_rmsd >= RMSthres:
                        out_mol2_mols.append(probe_mol)
                        out_mol2_blocks.append(probe_mol2_blocks[i])
                ref_mol2_mols = out_mol2_mols
            except:
                continue
    return out_mol2_blocks
