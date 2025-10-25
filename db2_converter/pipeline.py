import os
import subprocess
import shutil
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
from posebusters import PoseBusters
import logging

logger = logging.getLogger("DB2 generation")


from db2_converter.utils.rdkit_gen import rdk_enumerate_smi
from db2_converter.utils.utils import (
    next_mol2_lines,
    update_mol2block_from_mol,
    check_mol2_smi,
    run_external_command,
    exist_size,
    raise_errlog,
)
from db2_converter.utils.fixmol2 import fixmol2_by_template, fixmol2_and_du
from db2_converter.utils.prepare import prepare, create_namedata
from db2_converter.utils.conf_sample import (
    conf_sample,
    rdkit_prep,
    get_num_confs_for_mol,
)
from db2_converter.amsol.calc_charge_solvation import calc_charge_solvation
from db2_converter.utils.rmsd import RMSDfilter
from db2_converter.utils.match_frags import (
    f_AlignMolConformers,
    embed_blocks_molconformer,
    find_central_rigid,
    mol_to_ring_frags,
)
from db2_converter.mol2db2 import mol2db2
from db2_converter.mol2db2_py3_strain import mol2db2 as mol2db2_38
from db2_converter.strain.Mol2_Strain import mol2_strain
from db2_converter.config import config

# Basic config parameters
UNICON_EXE = config["all"]["UNICON_EXE"]
ANTECHAMBER = config["all"]["ANTECHAMBER"]

# pre-defined parameters
# Don't change if you do not know what they mean.
# EPS = 0.05  # RMSD likelihood threshold
EPS = 0.15
EPS_modify = 1e-4  # lower kept, upper modified
max_amsol_attempts = 5  # amsol attemp limit
RMSthres = 0.5  # RMSD threshold
sb_smarts = "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"
minfragsize = 3


def max_conf_assign(smi, max_conf, max_conf_noring):
    mol = Chem.MolFromSmiles(smi)
    if not len(mol_to_ring_frags(mol, cutsmarts=sb_smarts)):
        max_conf = max_conf_noring
    return max_conf


def sample_tp_unicon(infile,outfile=""):
    if not outfile:
        outfile = infile.stem + ".uni.smi"
        # only generate one smiles
        run_external_command(f"{UNICON_EXE} -i {infile} -o {outfile} -t single -p single")
    return Path(outfile)


def write_enumerated_smifile(inlines, outfile, method):
    logger.info(
        "####### Enumerating possible undefined stereochemistry of input SMILES... #######"
    )
    newlines = []
    toomany = False
    for line in inlines:
        try:
            smi, name = line.split()[0], line.split()[1]
            outsmis = rdk_enumerate_smi(smi)
            if len(outsmis) == 1:
                newlines.append(line)
            else:
                if len(outsmis) <= 2**5:
                    for i, outsmi in enumerate(outsmis):
                        newlines.append(f"{outsmi} {name}.{i}")
                        logger.info(f">>> {name}.{i} {outsmi}")
                else:
                    logger.error(
                        f">>> SMILES of {name} has too many stereoisomers, so skipped!"
                    )
                    toomany = True
                    faillist = [smi, name, method, "0enumerate_many"]
        except Exception as e:
            logger.error(e)
            logger.error(f">>> SMILES of {name} cannot be parsed")
    logger.info("####### Enumeration finished. #######\n")

    with open(outfile, "w") as f:
        for line in newlines:
            f.write(f"{line}\n")
    if toomany:
        return faillist


def fixmol2_wrapper(in_mol2lines_List, templateMol2lines="", smiles="", samplopt=""):
    tmp0mol2 = "tmp0.mol2"
    tmp0fixmol2 = "tmp0.mol2.fixed.mol2"
    TMPmol2 = "conformer.TMP.fixed.mol2"
    Path(tmp0mol2).write_text("".join(in_mol2lines_List[0]))
    # tentatively keeping antechamber converter
    run_external_command(
        f"{ANTECHAMBER} -i {tmp0mol2} -fi mol2 -o {tmp0fixmol2} -fo mol2 -at sybyl -pf y"
    )
    if not templateMol2lines:
        if not smiles:
            assert smiles != "", "SMILES is not given to fixmol2"
        if not exist_size(tmp0fixmol2) or not check_mol2_smi(tmp0fixmol2, smiles):
            # antechamber cannot deal with "[N-]" due to unexpected valence.
            shutil.copy(tmp0mol2, tmp0fixmol2)
        tmp0fix_mol2lines = open(tmp0fixmol2).readlines()
        tmp0fix_mol2lines_fix = fixmol2_and_du(smiles, tmp0fix_mol2lines)
        templateMol2lines = tmp0fix_mol2lines_fix
    if samplopt == "rdkit":
        # shutil.move(tmp0fixmol2, TMPmol2)
        Path(TMPmol2).write_text("".join(templateMol2lines))
    else:
        out_mol2lines_List = [
            fixmol2_by_template(mol2lines, templateMol2lines)
            for mol2lines in in_mol2lines_List
            ]
        return out_mol2lines_List


def chemistrycheck(insmi, in_mol2lines_List, checkstereo=True):
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*") # disable RDLogger of MolToInchi

    canonical_smiles = Chem.MolToSmiles(
        Chem.MolFromSmiles(insmi), isomericSmiles=checkstereo, canonical=True
    )
    all_mols = [
            Chem.MolFromMol2Block("".join(mol2lines), removeHs=True)
            for mol2lines in in_mol2lines_List
            ]
    out_mol2lines_List = []
    for i, mol in enumerate(all_mols):
        if checkstereo:
            # Chem.AssignAtomChiralTagsFromStructure(mol)
            """From rdkit:
            NOTE that this function does not check if atoms are chiral centers (i.e. all substituents are different),
            it merely sets the chiral type flags based on the coordinates and atom ordering.
            Use AssignStereochemistryFrom3D() if you want chiral flags only on actual stereocenters.
            """
            Chem.AssignStereochemistryFrom3D(mol)
        ref_inchi = Chem.MolToInchi(Chem.MolFromSmiles(canonical_smiles))
        gen_inchi = Chem.MolToInchi(mol)
        if gen_inchi == ref_inchi:
            out_mol2lines_List.append(in_mol2lines_List[i])
        else:
            logger.warning(f"Not consistent:\ngen {gen_inchi}\nref {ref_inchi}")
    logger.info(f">>> {len(out_mol2lines_List)} / {len(in_mol2lines_List)} passed chemistrycheck.")
    return out_mol2lines_List


def PB_filter(name, in_mol2lines_List):
    out_mol2lines_List = []
    PBmol2file = f"conformer.{name}.fixed.pb.mol2"
    with open(PBmol2file, "w") as f:
        for mol2lines in in_mol2lines_List:
            f.write("".join(mol2lines))
    run_external_command(
        f"{UNICON_EXE} -i {PBmol2file} -o {PBmol2file}.unipb.sdf",
    )
    buster = PoseBusters(config="mol")
    df = buster.bust([f"{PBmol2file}.unipb.sdf"], full_report=True)
    # df = buster.bust([f"{inmol2}"], full_report=True)
    df["PBvalid_conf"] = (
        df["all_atoms_connected"]
        & df["sanitization"]
        & df["bond_lengths"]
        & df["bond_angles"]
        & df["internal_steric_clash"]
        & df["aromatic_ring_flatness"]
        & df["double_bond_flatness"]
        & df["internal_energy"]
    )
    df = df.reset_index(drop=False)
    correct_indexes = df.index[df["PBvalid_conf"]].tolist()
    out_mol2lines_List = [in_mol2lines_List[i] for i in correct_indexes]
    logger.info(
        f">>> {len(out_mol2lines_List)} / {len(in_mol2lines_List)} passed PoseBusters filter."
    )
    return out_mol2lines_List


def mmffopt(in_mol2lines_List):
    out_mol2lines_List = []
    for mol2lines in in_mol2lines_List:
        newmol = Chem.MolFromMol2Block("".join(mol2lines), removeHs=False)
        MMFFOptimizeMolecule(newmol)
        out_mol2lines_List.append(update_mol2block_from_mol(mol2lines, newmol))
    return out_mol2lines_List


def prepare_mol2(name, smiles, inmol2, mergeiso=True):
    """
    # prapare.py
    ## generate a file called name.txt, one example content is:
    ## name.txt 1 CHEMBL85549 OC[C@H]1OC(N2C=NC3=C(NC4CCCCC4)N=CN=C32)[C@H](O)[C@@H]1O | NO_LONG_NAME
    ## This file will be used as default -n option in mol2db2.py
    """
    if mergeiso:
        name = name.split(".")[0]
    namedata = create_namedata(name=name, smiles=smiles, longname="NO_LONG_NAME")
    # for DUD-E / DUDE-Z when multiple protomers are considered as a single molecule
    # namedata = create_namedata(name=name.split(".")[0].split("_")[0], smiles=smiles, longname="NO_LONG_NAME")
    prepare(inmol2, namedata=namedata, namepath=None, tempdir=".")


def match_and_convert_mol2(
    mol2file,
    extra_fragsindex=[],
    extra_fragsmarts="",
    onlyextrafrags=False,
    reseth=True,
    rotateh=True,
    prefix="output",
    dock38=False,
):
    all_mol2lines = [x for x in next_mol2_lines(mol2file)]
    templateMol2lines = all_mol2lines[0]
    mol = Chem.MolFromMol2Block("".join(templateMol2lines), removeHs=False)
    # Get atom name to serial number mapping
    index_of_name = dict()
    for atom in mol.GetAtoms():
        name = atom.GetPropsAsDict()["_TriposAtomName"]  # extract O1, C2
        index_of_name[name] = atom.GetIdx()
    # Get break down fragments
    suitable_ring_frags = mol_to_ring_frags(mol=mol, cutsmarts=sb_smarts, minfragsize=3)
    fragsindex = []
    for frag in suitable_ring_frags:
        findex = []
        for atom in frag.GetAtoms():
            if atom.GetAtomicNum() > 1:  # Heavy atoms
                name = atom.GetPropsAsDict()["_TriposAtomName"]
                findex.append(index_of_name[name])
        # Double check to ensure at least 3 heavy atoms in one fragment
        if len(findex) >= minfragsize:
            fragsindex.append(findex)
    # extra rigid fragment #
    if extra_fragsindex:
        extra_fragsindex = [extra_fragsindex]
    if extra_fragsmarts:
        extra_fragsindex = []
        for extra_fragindex in mol.GetSubstructMatches(
            Chem.MolFromSmarts(extra_fragsmarts)
        ):
            extra_fragsindex += [list(extra_fragindex)]
    if extra_fragsindex:  # extra_fragsmarts has higher priority than extra_fragsindex
        fragsindex += extra_fragsindex
    if onlyextrafrags:
        fragsindex = extra_fragsindex
    ########################
    if not onlyextrafrags and not fragsindex: # means a molecule with no identified rigid part
        fragsindex = [find_central_rigid(mol)]
        logger.info(f">>> No rigid part can be kept after bond cleavage.")
        logger.info(
            f">>>>>> Will use central atom and its 1st neighbors {fragsindex}..."
        )
    if not fragsindex: # no rigid fragments can be found
        return -1
    mol_withconfs = embed_blocks_molconformer(all_mol2lines, removeHs=False)
    # Aligning every frag
    i = 0
    for dirname in ["sdf", "mol2", "db2"]:
        Path(dirname).mkdir(exist_ok=True)
    for index in fragsindex:
        all_index = f_AlignMolConformers(
            mol_withconfs,
            atomIds=index,
            maxIters=500,
            reflect=False,
            eps=EPS,
            eps_modify=EPS_modify,
        )  # maxIters and reflect are for subfunc AlignMolConformers
        # all_index is a list of lists, each list with len>1 means a cluster
        logger.debug("Rigid body index: %s", index)
        logger.debug("Conformers used id: %s", all_index)
        for conf_index in all_index:
            writer = Chem.SDWriter(f"sdf/{prefix}.{i}.sdf")
            for conf in conf_index:
                writer.write(
                    mol_withconfs, confId=conf
                )  # write different aligned confs into separated sdf files
            writer.close()
            run_external_command(
                f"{UNICON_EXE} -i sdf/{prefix}.{i}.sdf -o  mol2/{prefix}.{i}.mol2",
                stderr=subprocess.STDOUT,
            )
            # if output mol2 has format issue, put fixmol2_wrapper here
            in_mol2lines_List = list(next_mol2_lines(f"mol2/{prefix}.{i}.mol2"))
            out_mol2lines_List = fixmol2_wrapper(in_mol2lines_List, templateMol2lines)
            with open(f"mol2/{prefix}.{i}.fixed.mol2", "w") as f:
                for mol2lines in out_mol2lines_List:
                    f.write("".join(mol2lines))
            if not dock38:
                ############## DOCK37 db2 from mol2 ##############
                mol2db2.mol2db2_main(
                    mol2file=f"mol2/{prefix}.{i}.fixed.mol2",
                    solvfile=f"{prefix}.solv",
                    namefile="name.txt",
                    db2gzfile=f"db2/{prefix}.{i}.db2.gz",
                    timeit=True,
                    reseth=reseth,
                    rotateh=rotateh,
                )
                ##################################################
            else:
                ############## DOCK38 db2 from mol2 ##############
                db2in_standard = mol2_strain(f"mol2/{prefix}.{i}.fixed.mol2")
                mol2db2_38.mol2db2_main(
                    db2in_standard,
                    solvfile=f"{prefix}.solv",
                    db2gzfile=f"db2/{prefix}.{i}.db2.gz",
                    timeit=True,
                    reseth=reseth,
                    rotateh=rotateh,
                )
            i += 1
    return i


def gen_conf(
    zinc,
    max_conf,
    max_conf_noring,
    samplopt,
    checkstereo,
    PBfilter=False,
    MMFFopt=False,
    cluster=False,
    limitconf=False,
    extra_fragsindex=[],
    extra_fragsmarts="",
    onlyextrafrags=False,
    keep_max_conf=False,
    reseth=True,
    rotateh=True,
    prefix="output",
    faillist=[],
    bcl_option="",
    confgenx_option="",
    mergeiso=True,
    dock38=False,
    **kwargs
):
    logger.info(f"############### Now dealing with {zinc}... ###############")
    Path(zinc).mkdir(exist_ok=True)
    shutil.copy(f"{zinc}.smi", zinc)
    os.chdir(zinc)
    smi = Path(f"{zinc}.smi").read_text().split("\n")[0].split()[0]
    ringmol = Chem.MolFromSmiles(smi)
    N_ring = len(mol_to_ring_frags(ringmol, cutsmarts=sb_smarts))
    if rotateh:
        NrotHs, multiplier = mol2db2.mol2db2_to_numhyds(f"{zinc}.smi")
    if not keep_max_conf:
        if not N_ring and not (extra_fragsmarts or extra_fragsindex):
            max_conf = max_conf_noring
        if N_ring and rotateh:  # if rotateh, we need to truncate the ensemble size
            if NrotHs in [4, 5]:
                max_conf = max_conf // 30
            if NrotHs in [2, 3]:
                max_conf = max_conf // 3
        if limitconf:
            max_conf = min(get_num_confs_for_mol(smi), max_conf)
    elif limitconf:
        max_conf = min(get_num_confs_for_mol(smi), max_conf)
    # if keep_max_conf, the max_conf you defined is the actual max_conf for generation
    N_max_conf_in = max_conf

    samplname = samplopt
    samplopts = samplopt.split("_")
    if len(samplopts) > 1:
        logger.info(f"#### {samplname} ensemble sampling... ####")

    mol2lines_List_dict = {}
    for samplopt in samplopts:
        mol2file = f"conformer.{zinc}.{samplopt}.mol2"
        fixed_mol2file = f"conformer.{zinc}.{samplopt}.fixed.mol2"

        logger.info(
            f">>> {samplopt} used to sample {max_conf} conformers, checkstereo {checkstereo}, PBfilter {PBfilter}, MMFFopt {MMFFopt}, cluster {cluster}"
        )
        if samplopt == "rdkit":
            rdkit_prep(zinc, mol2file)
        else:
            conf_sample(
                samplopt,
                zinc,
                mol2file,
                max_conf,
                bcl_option,
                confgenx_option,
                log=logger,
            )

        if not exist_size(mol2file):
            error = "1generate"
            faillist.append([smi, zinc, samplopt, error])
            raise_errlog(error, logger)
            mol2lines_List_dict[samplopt] = []
            continue

        try:
            in_mol2lines_List = list(next_mol2_lines(mol2file))
            out_mol2lines_List = fixmol2_wrapper(in_mol2lines_List, "", smi, samplopt)
        except Exception as e:
            logger.error(e)
            error = "2fixmol2"
            faillist.append([smi, zinc, samplopt, error])
            raise_errlog(error, logger)
            mol2lines_List_dict[samplopt] = []
            continue

        if samplopt == "rdkit":
            try:
                conf_sample(samplopt, zinc, fixed_mol2file, max_conf, log=logger)
                out_mol2lines_List = list(next_mol2_lines(fixed_mol2file))
            except Exception as e:
                logger.error(e)
                error = "1generate"
                faillist.append([smi, zinc, samplopt, error])
                raise_errlog(error, logger)
                mol2lines_List_dict[samplopt] = []
                continue

        if not out_mol2lines_List:
            error = "2fixmol2"
            faillist.append([smi, zinc, samplopt, error])
            raise_errlog(error, logger)
            mol2lines_List_dict[samplopt] = []
            continue
        in_mol2lines_List = out_mol2lines_List

        ## Chemistry check and filter (Through RDKit)
        if checkstereo:
            logger.info(">>> Chemistry checking with stereochemistry...")
        else:
            logger.info(">>> Chemistry checking without stereochemistry... Be careful!!!")
        try:
            out_mol2lines_List = chemistrycheck(
                insmi=smi,
                in_mol2lines_List=in_mol2lines_List,
                checkstereo=checkstereo,
            )
        except Exception as e:
            out_mol2lines_List = []
        if not out_mol2lines_List:
            error = "3chemistrycheck"
            faillist.append([smi, zinc, samplopt, error])
            raise_errlog(error, logger)
            mol2lines_List_dict[samplopt] = []
            continue
        else:
            logger.info(">>> Chemistry Checked.")
        in_mol2lines_List = out_mol2lines_List

        ## PoseBusters filter
        if PBfilter:
            logger.info(">>> PoseBusters filter implausible conformers...")
            out_mol2lines_List = PB_filter(name=zinc, in_mol2lines_List=in_mol2lines_List)
            logger.info(">>> PoseBusters filter Finished.")
        if not out_mol2lines_List:
            error = "4PBfilter"
            faillist.append([smi, zinc, samplopt, error])
            raise_errlog(error, logger)
            mol2lines_List_dict[samplopt] = []
            continue
        in_mol2lines_List = out_mol2lines_List

        ## MMFF optimization
        if MMFFopt:
            logger.info(">>> MMFF optimization...")
            out_mol2lines_List = mmffopt(in_mol2lines_List)
            logger.info(">>> MMFF optimization Finished.")

        ## Iterate to get one normal amsol output
        logger.info(">>> AMSOL calculation for atom partial charges and desolvation...")
        if not exist_size(f"{prefix}.solv"):  # keep for ensemble
            for i, mol2lines in enumerate(out_mol2lines_List):
                logger.debug(f"amsol trial {i}")
                with open(f"{prefix}{zinc}.mol2", "w") as f:
                    f.write("".join(mol2lines))

                prepare_mol2(
                    name=zinc,
                    smiles=smi,
                    inmol2=f"{prefix}{zinc}.mol2",
                    mergeiso=mergeiso,
                )

                # AMSOL ( very complicated in it, but not much work )
                ## with a lot of output log info
                ## ::input: output$number.mol2 (single mol2), total charge used by AMSOL is directly summed read from the APC in it
                ## ::output: {prefix}.solv (amsol output file)
                ## ::output: {prefix}.mol2 (single mol2, w/ charge)
                ## ::output: {prefix}$number.mol2 (single mol2, format slightly changed, w/o charge)
                try:
                    calc_charge_solvation(f"{prefix}{zinc}.mol2")
                except:
                    pass
                if exist_size(f"{prefix}.solv") or i + 1 == max_amsol_attempts:
                    logger.info(">>> AMSOL calculation Finished.")
                    break

            if not exist_size(f"{prefix}.solv"):
                error = "5amsolfail"
                faillist.append([smi, zinc, samplopt, error])
                raise_errlog(error, logger)
                continue
        faillist.append([])
        # If success, faillist will be [], this could be helpful for combinatorial generation

        mol2lines_List_dict[samplopt] = out_mol2lines_List

    # collect mol2, rmsd filter
    in_mol2lines_List = []
    fixed_mol2file = f"conformer.{zinc}.fixed.mol2"
    fail_name_count = len([onefail for onefail in faillist if onefail])
    if fail_name_count == len(samplopts):
        return faillist, N_max_conf_in, "", "", N_ring, ""
    else:
        for samplopt in samplopts:
            in_mol2lines_List += mol2lines_List_dict[samplopt]
        if not cluster:
            out_mol2lines_List = in_mol2lines_List
        else:
            logger.info(f">>> RMSD-based Clustering at {RMSthres} Angstrom...")
            out_mol2lines_List = RMSDfilter(mol2lines_List_dict, zinc, samplopts, RMSthres)
            logger.info(f">>> {len(out_mol2lines_List)} / {len(in_mol2lines_List)} kept.")
            logger.info(">>> RMSD Clustering finished...")
        with open(fixed_mol2file, "w") as f:
            f.write("".join(["".join(mol2_block) for mol2_block in out_mol2lines_List]))
    N_act_conf = len(out_mol2lines_List)
    if rotateh:
        N_act_conf_out_rotH = N_act_conf * multiplier
    else:
        N_act_conf_out_rotH = N_act_conf

    # Match and convert
    ## ::input: conformer.$number.fixed.mol2 (infile)(multiple conformations)
    ## ::output: dir sdf with separated cluster sdf file
    ## ::output: dir mol with separated cluster mol2 file
    ## ::output: dir db2 with separated db2.gz file
    shutil.rmtree("db2", ignore_errors=True)
    db2part_count = match_and_convert_mol2(
        mol2file=fixed_mol2file,
        extra_fragsindex=extra_fragsindex,
        extra_fragsmarts=extra_fragsmarts,
        onlyextrafrags=onlyextrafrags,
        reseth=reseth,
        rotateh=rotateh,
        prefix="output",
        dock38=dock38
    )
    # collect
    if Path("db2").exists():
        subprocess.run("cat db2/*.db2.gz > all.db2.gz", shell=True)
    if os.path.getsize("all.db2.gz") != 28:
        logger.info(f">>> {db2part_count} db2 blocks were generated and saved...")
        logger.info(f"############### Finished with {zinc}... ###############\n")
    else:
        error = "9nulldb2gz"
        faillist.append([smi, zinc, samplopt, error])
        raise_errlog(error, logger)

    return (
        faillist,
        N_max_conf_in,
        N_act_conf,
        N_act_conf_out_rotH,
        N_ring,
        db2part_count,
    )
