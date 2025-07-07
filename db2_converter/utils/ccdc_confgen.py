def activate_ccdc():
    import logging
    import sys
    import subprocess
    from db2_converter.config import config
    from db2_converter.utils.utils import run_external_command
    logger = logging.getLogger("db2_converter")
    try:
        CCDC_PYTHON3 = config["ccdc"]["CCDC_PYTHON3"]
        CCDC_activate_node = config["ccdc"]["CCDC_activate_node"]
        CCDC_activate_command = config["ccdc"]["CCDC_activate_command"]
        CCDC_activate_code = config["ccdc"]["CCDC_activate_code"]
        ccdc_activate_whole_command = f"{CCDC_activate_command} -a -k {CCDC_activate_code}"
        ccdc_license_check_command = f"{CCDC_PYTHON3} -c 'import ccdc'"
        if not CCDC_activate_node in [ "local", "localhost" ]:
            ccdc_activate_whole_command = f"ssh {CCDC_activate_node} {ccdc_activate_whole_command}"
            run_external_command(ccdc_activate_whole_command)
        else:
            if run_external_command(ccdc_license_check_command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE): # if no output, already activated
                run_external_command(ccdc_activate_whole_command) # otherwise, should re-activate
    except Exception:
        logger.error("!!! CCDC is called but activation has failed...")
        sys.exit()


if __name__ == "__main__":
    from ccdc import conformer
    from ccdc import io
    import argparse
    parser = argparse.ArgumentParser("python API of CSD Conformer Generator")
    parser.add_argument("--infile", type=str, required=True, help="input molecule file, support cif|csdsql|csdsqlx|identifiers|mol|mol2|res|sdf|sqlite|sqlmol2, should be protonated")
    parser.add_argument("--outfile", type=str, required=True, help="output molecule file, support cif|csdsql|identifiers|mol|mol2|pdb|res|sdf")
    parser.add_argument("--max_conf", type=int, required=True, help="maximum number of conformers that can be genearted")
    parser.add_argument("--nthreads", type=int, default=4, help="number of cpu cores used")
    parser.add_argument("--max_unusual_torsions", type=int, default=2, help="Conformer Generator Setting:::max_unusual_torsions")
    parser.add_argument('--lock_bond_list', nargs='+', default=[], help="Should be 2x list, for example C13 C14, or C11, C12, C13, C14")
    args = parser.parse_args()

    infile = args.infile
    outfile = args.outfile
    max_conf = args.max_conf
    nthreads = args.nthreads
    max_unusual_torsions = args.max_unusual_torsions 
    lock_bond_list = args.lock_bond_list

    mol_reader = io.MoleculeReader(infile)
    mol = mol_reader[0] # ccdc.molecule.Molecule
    
    conformer_generator = conformer.ConformerGenerator(
        settings=None,
        skip_minimisation=False, # Initialize the input 3D molecular structure's bond lengths and valence angles
        nthreads=nthreads,
        parameter_locator=conformer.DefaultConformerParameterFileLocator()
    )

    # Default values
    # print(conformer_generator.settings.max_conformers) # 200
    # print(conformer_generator.settings.max_unusual_torsions) # 2, meaning conformations with more than 2 unusual rotamers are disallowed.
    # print(conformer_generator.settings.superimpose_conformers_onto_reference) # False
    # print(conformer_generator.settings._use_input_torsion_distributions) # False
    """ Conformer Generator Settings """
    conformer_generator.settings.max_conformers = max_conf 
    conformer_generator.settings.max_unusual_torsions = 2
    conformer_generator.settings.superimpose_conformers_onto_reference = False
    conformer_generator.settings._use_input_torsion_distributions = False

    for i in range(0, len(lock_bond_list), 2):
        lock_bond_atom1, lock_bond_atom2 = lock_bond_list[i], lock_bond_list[i+1]
        conformer_generator.lock_torsion(mol.bond(lock_bond_atom1, lock_bond_atom2)) # Optional, lock torsion bond
        # Specify that a particular torsion should not be changed when generating conformers of its molecule.
        # If the bond is in a ring, the whole ring will be locked.
    
    """ Conformer Generator """
    conformers = conformer_generator.generate(mol)
    print(f">>> {len(conformers)} / {max_conf} conformers generated.")

    print(">>>>>> Whether the internal sampling limit as been reached:", conformers.sampling_limit_reached)
    print(">>>>>> Number of rotamers for which no crystallographic data is available.", conformers.n_rotamers_with_no_observations)
    
    """ Write Conformer Ensemble into a Molecule """
    conformers_mols = [ c.molecule for c in conformers ]
    with io.MoleculeWriter(outfile) as mol_writer:
            for mol in conformers_mols:
                    mol_writer.write(mol)
