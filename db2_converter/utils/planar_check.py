# posebusters.modules.flatness
from copy import deepcopy
from typing import Any, Iterable

import numpy as np
from numpy import ndarray as Array
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolTransforms import GetDihedralRad
from rdkit.Chem.rdmolfiles import MolFromSmarts

_flat = {
    "aromatic_5_membered_rings_sp2": "[ar5^2]1[ar5^2][ar5^2][ar5^2][ar5^2]1",
    "aromatic_6_membered_rings_sp2": "[ar6^2]1[ar6^2][ar6^2][ar6^2][ar6^2][ar6^2]1",
    "trigonal_planar_double_bonds": "[C;X3;^2](*)(*)=[C;X3;^2](*)(*)",
    "amide_bonds": "[CX3](=[OX1])[NX3][H]"
}


def check_flatness(
    mol_pred: Mol, threshold_flatness: float = 0.1, flat_systems: dict[str, str] = _flat
) -> dict[str, Any]:
    """Check whether substructures of molecule are flat.

    Args:
        mol_pred: Molecule with exactly one conformer.
        threshold_flatness: Maximum distance from shared plane used as cutoff. Defaults to 0.1.
        flat_systems: Patterns of flat systems provided as SMARTS. Defaults to 5 and 6 membered
            aromatic rings and carbon sigma bonds.

    Returns:
        PoseBusters results dictionary.
    """
    mol = deepcopy(mol_pred)

    planar_groups = []
    types = []
    for flat_system, smarts in flat_systems.items():
        match = MolFromSmarts(smarts)
        atom_groups = list(mol.GetSubstructMatches(match))
        if flat_system == "amide_bonds":
            for atom_group in atom_groups:
                Cidx, Oidx, Nidx, Hidx = atom_group
                conf = mol.GetConformer()
                if np.cos(GetDihedralRad(conf,Oidx,Cidx,Nidx,Hidx)) > 0:
                    return False
        planar_groups += atom_groups
        types += [flat_system] * len(atom_groups)

    # calculate distances to plane and check threshold
    coords = [_get_coords(mol, group) for group in planar_groups]
    max_distances = [float(_get_distances_to_plane(X).max()) for X in coords]
    flatness_passes = [bool(d <= threshold_flatness) for d in max_distances]

    result = all(flatness_passes) if len(flatness_passes) > 0 else True

    return result


def _get_distances_to_plane(X: Array) -> Array:
    """Get distances of points X to their common plane."""
    # center points X in R^(n x 3)
    X = X - X.mean(axis=0)
    # singular value decomposition
    _, _, V = np.linalg.svd(X)
    # last vector in V is normal vector to plane
    n = V[-1]
    # distances to plane are projections onto normal
    d = np.dot(X, n)
    return d


def _get_coords(mol: Mol, indices: Iterable[int]) -> Array:
    return np.array([mol.GetConformer().GetAtomPosition(i) for i in indices])
