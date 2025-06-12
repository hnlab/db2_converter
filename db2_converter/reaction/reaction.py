import sys
import logging
import xml.etree.ElementTree as ET
from collections import namedtuple
from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger("Reaction")

REAOTHERCOLOR=21
CAPCOLOR=22
DEFAULTCOLOR=23
LINKHEAVYCOLOR=24
LINKHYDROCOLOR=25

ReactionInfoClass = namedtuple(
    "ReactionInfo",
    [
        "ncore",
        "ncap",
        "n_reagents",
        "n_products",
        "reagent2_smi",
        "reactionSMARTS",
        "productSMARTS",
        "protectidx_list",
        "protectSMARTS_list",
        "filterSMARTS_list",
        "chemcolor",
        "neutralizeAtomId_list"
    ]
)

def parse_reaction_xmlfile(xmlfile,tgt_reaction_name,tgt_reagent_type,tgt_cap):
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    reaction_matched = False
    for reactionClass in root.findall("reactionClass"):
        reactionClass_name = reactionClass.get("name")
        if not reactionClass_name == tgt_reaction_name:
            continue
        reagentTypeList = reactionClass.find("reagentTypeList")
        reagenttypes = [ reagent.get("name") for reagent in reagentTypeList.findall(f"reagent") ]
        if not tgt_reagent_type in reagenttypes:
            continue
        reactionList = reactionClass.find("reactionList")
        reactions = reactionList.findall("reaction")
        for reaction in reactions:
            this_cap = int(reaction.get("cap"))
            this_reagenttype = reaction.get("reagenttype")
            if tgt_cap == this_cap and tgt_reagent_type == this_reagenttype:
                reaction_matched = True
                reactionSMARTS = reaction.get("reactionSMARTS")
                reagent2_smi = reaction.get("reagent2smi")
                reaction_matched = True
                break
    if reaction_matched:
        # reactionInfo
        reactionInfo = reactionClass.find("reactionInfo")
        core = int(reactionInfo.get("core"))
        productSMARTS = reactionInfo.get("productSMARTS")
        n_reagents = int(reactionInfo.get("reagentcount"))
        n_products = int(reactionInfo.get("productcount"))
        chemcolor = [ int(i) for i in reactionInfo.get("chemcolor").split() ]
        # reactionInfo
        # neutralizeList
        neutralizeList = reactionClass.find("neutralizeList")
        neutralizeAtomId_list = [
             int(atom.get("index"))
             for atom in neutralizeList.findall("atom")
             if atom.get("reagenttype") == tgt_reagent_type
        ]
        # neutralizeList
        # protect
        protectList = reactionClass.find("protectList")
        protects = [ protect for protect in protectList.findall("protect") if protect.get("reagenttype") == tgt_reagent_type]
        protectidx_list = [ int(protect.get("protectidx")) for protect in protects ]
        protectSMARTS_list = [ protect.get("protectSMARTS") for protect in protects ]
        # protect
        # filter
        filterList = reactionClass.find("filterList")
        filterSMARTS_list = [ onefilter.get("filterSMARTS") for onefilter in filterList.findall("filter") ]
        # filter
        reaction_info_list = [
            core,
            tgt_cap,
            n_reagents,
            n_products,
            reagent2_smi,
            reactionSMARTS,
            productSMARTS,
            protectidx_list,
            protectSMARTS_list,
            filterSMARTS_list,
            chemcolor,
            neutralizeAtomId_list
            ]
        reactInfoObj = ReactionInfoClass(*reaction_info_list)
        return reactInfoObj
    else:
        logger.error("Reaction not found!!!\nCarefully check your reaction input parameter!!!\nGeneration failed and exit.")
        # reaction error log
        sys.exit()


def reaction(reagent1_smi,reactInfoObj):
    rxn = AllChem.ReactionFromSmarts(reactInfoObj.reactionSMARTS)
    reagent1_mol = Chem.MolFromSmiles(reagent1_smi)
    reagent2_mol = Chem.MolFromSmiles(reactInfoObj.reagent2_smi)
    # set protection
    for protectidx, protectSMARTS in zip(reactInfoObj.protectidx_list, reactInfoObj.protectSMARTS_list):
        for match in reagent1_mol.GetSubstructMatches(Chem.MolFromSmarts(protectSMARTS)):
            reagent1_mol.GetAtomWithIdx(match[protectidx]).SetProp('_protected','1')
    product = rxn.RunReactants((reagent1_mol,reagent2_mol))
    # reaction
    for oneproduct in product:
        prodct_mol = oneproduct[0] # only 1 product
        prodct_smi = Chem.MolToSmiles(prodct_mol)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(prodct_mol)) # only 1 product
        react_idxs = list(range(reactInfoObj.ncore)) # react idxs, e.g., 0,1,2,3
        cap_idxs = list(range(reactInfoObj.ncore, reactInfoObj.ncore + reactInfoObj.ncap)) # capped methyl C e.g., 4
        cano_map_index = map_index(prodct_mol)
        react_idxs = [ cano_map_index[idx] for idx in react_idxs ]
        cap_idxs = [ cano_map_index[idx] for idx in cap_idxs ]
        # neutralize
        neutralizeAtomIds = [ cano_map_index[neutralizeAtomId] for neutralizeAtomId in reactInfoObj.neutralizeAtomId_list ]
        mol = neutralize_atoms(mol,neutralizeAtomIds)
        # neutralize
        if not mol:
            continue # neutralize failed.
        if not mol.GetSubstructMatches(Chem.MolFromSmarts(reactInfoObj.productSMARTS)): # check if generate the expected product
            continue
        # check
        for filterSMARTS in reactInfoObj.filterSMARTS_list:
            if mol.GetSubstructMatches(Chem.MolFromSmarts(filterSMARTS)):
                mol = None
                break
        if mol:
            smi = Chem.MolToSmiles(mol)
            chem_color_dict = dict(zip(react_idxs,reactInfoObj.chemcolor))
            for idx in cap_idxs:
                chem_color_dict[idx] = CAPCOLOR
            yield (smi,react_idxs,chem_color_dict)


def neutralize_atoms(mol,atomIdxs):
    # source: https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules
    for at_idx in atomIdxs:
        atom = mol.GetAtomWithIdx(at_idx)
        chg = atom.GetFormalCharge()
        if chg != 0: # only operate atoms with formal charge
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            if hcount - chg < 0:
                return None
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def map_index(mol):
    return dict(
        zip(
            list(map(int, mol.GetProp("_smilesAtomOutputOrder")[1:-2].split(","))),
            list(range(mol.GetNumAtoms()))
        )
    )
