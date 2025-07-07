import sys
from pathlib import Path
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

XMLFILE = Path(__file__).parent / "reaction.xml"
def parse_reaction_xmlfile(xmlfile,tgt_reaction_name,tgt_reagent_type,tgt_cap):
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    reaction_matched = False
    for reactionClass in root.findall("reactionClass"):
        reactionClass_name = reactionClass.get("name")
        if not reactionClass_name == tgt_reaction_name:
            continue
        reactionInfo = reactionClass.find("reactionInfo")
        reagentTypeList = reactionClass.find("reagentTypeList")
        reagenttypes = [ reagent.get("name") for reagent in reagentTypeList.findall(f"reagent") ]
        if not tgt_reagent_type in reagenttypes:
            continue
        reactionList = reactionClass.find("reactionList")
        reactions = reactionList.findall("reaction")
        for reaction in reactions:
            this_cap = int(reaction.get("cap"))
            this_reagenttype = reaction.get("reagenttype")
            # print(tgt_cap == this_cap, tgt_reagent_type == this_reagenttype)
            if tgt_cap == this_cap and tgt_reagent_type == this_reagenttype:
                reaction_matched = True
                reactionSMARTS = reaction.get("reactionSMARTS")
                reagent2_smi = reaction.get("reagent2smi")
                break
        if reaction_matched:
            break
    if reaction_matched:
        # reactionInfo
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
        last_smi = ""
        ncore_plus_ncap = reactInfoObj.ncore + reactInfoObj.ncap
        prodct_mol = oneproduct[0] # only 1 product

        prodct_smi = Chem.MolToSmiles(prodct_mol, rootedAtAtom=ncore_plus_ncap-1)
        logger.debug(f"prodct_smi: {prodct_smi}")
        mol = Chem.MolFromSmiles(prodct_smi) # kekulize it
        if not mol:
            continue
        react_idxs = list(range(reactInfoObj.ncore)) # react idxs, e.g., 0,1,2,3
        cap_idxs = list(range(reactInfoObj.ncore, ncore_plus_ncap)) # capped methyl C e.g., 4
        logger.debug(f"react_idxs: {react_idxs}")
        logger.debug(f"cap_idxs: {cap_idxs}")
        cano_map_index = map_index(prodct_mol)
        react_idxs = [ cano_map_index[idx] for idx in react_idxs ]
        cap_idxs = [ cano_map_index[idx] for idx in cap_idxs ]
        logger.debug(f"map_react_idxs: {react_idxs}")
        logger.debug(f"map_cap_idxs: {cap_idxs}")
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
            smi = Chem.MolToSmiles(mol, rootedAtAtom=0) # keep the input orientation
            chem_color_dict = dict(zip(react_idxs,reactInfoObj.chemcolor))
            for idx in cap_idxs:
                chem_color_dict[idx] = CAPCOLOR
            if last_smi != smi: # remove multiple redundant matches
                logger.debug(f"output_prodct_smi: {smi}")
                yield (smi,react_idxs,chem_color_dict)
                last_smi = smi


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


def parse_reverse_reaction_dict(xmlfile):
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    reaction_group_dict = {}
    # "amide": "[CX3](=[OX1])[NX3]",
    # "ester": "[CX3](=[OX1])[OX2]",
    # "urea": "[NX3][CX3](=[OX1])[NX3]",
    # "thiourea": "[NX3][CX3](=[SX1])[NX3]",
    # "sulfoamide": "[SX4](=[OX1])(=[OX1])[NX3]",
    reagent_dict = {}
    # "amide": ["carboxylate", "amine"],
    # "ester": ["carboxylate", "alcohol"],
    # "urea": ["isocyanate", "amine"],
    # "thiourea": ["isothiocyanate", "amine"],
    # "sulfoamide": ["sulfonyl_chloride", "amine"],
    reverse_reaction_dict = {}
    # "amide": "[!$(*#*)&!D1:0][CX3:1](=[OX1:2])[NX3:3][!$(*#*)&!D1:8]>>[!$(*#*)&!D1:0][CX3:1](=[OX1:2])[OX2].[NX3:3][!$(*#*)&!D1:8]",
    # "ester": "[!$(*#*)&!D1:0][CX3:1](=[OX1:2])[OX2:3][!$(*#*)&!D1:8]>>[!$(*#*)&!D1:0][CX3:1](=[OX1:2])[OX2].[OX2:3][!$(*#*)&!D1:8]",
    # "urea": "[!$(*#*)&!D1:0][NX3:1][CX3:2](=[OX1:3])[NX3:4][!$(*#*)&!D1:8]>>[!$(*#*)&!D1:0][NX2:1]=[CX2:2]=[OX1:3].[NX3:4][!$(*#*)&!D1:8]",
    # "thiourea": "[!$(*#*)&!D1:0][NX3:1][CX3:2](=[SX1:3])[NX3:4][!$(*#*)&!D1:8]>>[!$(*#*)&!D1:0][NX2:1]=[CX2:2]=[SX1:3].[NX3:4][!$(*#*)&!D1:8]",
    # "sulfoamide": "[!$(*#*)&!D1:0][SX4:1](=[OX1:2])(=[OX1:3])[NX3:4][!$(*#*)&!D1:8]>>[!$(*#*)&!D1:0][SX4:1](=[OX1:2])(=[OX1:3])[Cl].[NX3:4][!$(*#*)&!D1:8]",


    for reactionClass in root.findall("reactionClass"):
        reactionName = reactionClass.get("name")
        reactionInfo = reactionClass.find("reactionInfo")
        productSMARTS = reactionInfo.get("productSMARTS")
        reagentTypeList = reactionClass.find("reagentTypeList")
        reagentNames = [ reagent.get("name") for reagent in reagentTypeList.findall("reagent") ]
        reverseReaction = reactionClass.find("reverseReaction")
        reverseReactionSMARTS = reverseReaction.get("reactionSMARTS")
        reverseProductSMARTS = reverseReactionSMARTS.split(">>")[0]
        reverseReagentSMARTS_list = reverseReactionSMARTS.split(">>")[1].split(".")
        nonterm_reverseReactionSMARTS = f"[!$(*#*)&!D1:0]{reverseProductSMARTS}[!$(*#*)&!D1:8]>>[!$(*#*)&!D1:0]{reverseReagentSMARTS_list[0]}.{reverseReagentSMARTS_list[1]}[!$(*#*)&!D1:8]"

        reaction_group_dict[reactionName] = productSMARTS
        reagent_dict[reactionName] = reagentNames
        reverse_reaction_dict[reactionName] = nonterm_reverseReactionSMARTS

    return reaction_group_dict, reagent_dict, reverse_reaction_dict
