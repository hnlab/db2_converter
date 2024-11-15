#!/usr/bin/env python
import sys
import logging

logger = logging.getLogger("amsol")


def read_ZmatMOPAC(Zmat_file):

    logger.debug("just entered read_ZmatMOPAC()")

    infile_ZmatMOPAC = open(Zmat_file, "r")
    lines = infile_ZmatMOPAC.readlines()

    ZmatMOPAC_lines = {}
    line_key_infile = 0
    line_key_out = 0
    spl = []

    # loop over all the lines in infile_ZmatMOPAC, whose first 3 lines
    # contain information which can be disregarded
    for line in lines:
        line_key_infile += 1
        # cut the first three lines off. They must be disregarded.
        if not (line_key_infile < 4):
            line_key_out += 1  # after cutting the first 3 lines, the fourth line should be the new first line and so forth.
            #
            # in the first three lines of the Z-matrix, obabel writes 1 too often. This would rise a non-fatal error message in amsol7.1.
            # this error message should be avoided, therefore the first three Z-matrix lines should be slightly modified
            if line_key_out == 1:
                spl = line.split()
                spl[2] = "0"  # instead of 1
                spl[4] = "0"  # instead of 1
                spl[6] = "0"  # instead of 1
                line = "%-2s %10.6f %2d %11.6f %2d %11.6f %2d %5d %3d %3d\n" % (
                    spl[0],
                    float(spl[1]),
                    int(spl[2]),
                    float(spl[3]),
                    int(spl[4]),
                    float(spl[5]),
                    int(spl[6]),
                    int(spl[7]),
                    int(spl[8]),
                    int(spl[9]),
                )
            if line_key_out == 2:
                spl = line.split()
                spl[4] = "0"  # instead of 1
                spl[6] = "0"  # instead of 1
                line = "%-2s %10.6f %2d %11.6f %2d %11.6f %2d %5d %3d %3d\n" % (
                    spl[0],
                    float(spl[1]),
                    int(spl[2]),
                    float(spl[3]),
                    int(spl[4]),
                    float(spl[5]),
                    int(spl[6]),
                    int(spl[7]),
                    int(spl[8]),
                    int(spl[9]),
                )
            if line_key_out == 3:
                spl = line.split()
                spl[6] = "0"  # instead of 1
                line = "%-2s %10.6f %2d %11.6f %2d %11.6f %2d %5d %3d %3d\n" % (
                    spl[0],
                    float(spl[1]),
                    int(spl[2]),
                    float(spl[3]),
                    int(spl[4]),
                    float(spl[5]),
                    int(spl[6]),
                    int(spl[7]),
                    int(spl[8]),
                    int(spl[9]),
                )

            ZmatMOPAC_lines[line_key_out] = line
    infile_ZmatMOPAC.close()

    logger.debug("read_ZmatMOPAC() has finished.")

    return ZmatMOPAC_lines


def create_amsol71_inputfile(Path_and_NameZmatMOPACFile, MoleculeName, ZmatMOPAC_Data):

    logger.debug("just entered the function create_amsol71_inputfile(): ")
    logger.debug("in the case of problems:")
    logger.debug("Make sure that OpenEye OEChem is installed on your system.")

    string_Path_And_NameZmatMOPACFile = Path_and_NameZmatMOPACFile

    # slice off the ending .ZmatMOPAC from string_Path_And_NameZmatMOPACFile by using [0:-10]:

    # open a file for the SM5.42R calculation in water solvent:
    Actual_Amsol71_InputFile_Water = open(
        "%s" % (string_Path_And_NameZmatMOPACFile[0:-10] + ".in-wat"), "w"
    )

    # open a file for the SM5.42R calculation in hexadecane solvent:
    Actual_Amsol71_InputFile_Hexadecane = open(
        "%s" % (string_Path_And_NameZmatMOPACFile[0:-10] + ".in-hex"), "w"
    )

    ## Modified by qiuyu
    infile = string_Path_And_NameZmatMOPACFile[0:-10] + ".mol2"
    atomlines = list()
    flag = False
    for line in open(infile):
        if "@<TRIPOS>ATOM" in line:
            flag = True
            continue
        elif "@<TRIPOS>BOND" in line:
            flag = False
            continue
        if flag:
            if line.strip():
                atomlines.append(line)
    netcharge = 0.0
    for line in atomlines:
        netcharge += float(line.split()[8])

    netcharge = 1.0 * round(netcharge)
    logger.debug(
        f"netcharge of molecule in temp.mol2 (sum of partial charges): {netcharge}"
    )

    # write the AMSOL7.1 keywords for a SM5.42R point calculation in water to the AMSOL7.1 water input-file:
    Water_Amsol71_SM542R_Keywords = (
        """CHARGE=%s AM1 1SCF TLIMIT=15 GEO-OK SM5.42R\n& SOLVNT=WATER\n""" % netcharge
    )
    Actual_Amsol71_InputFile_Water.write(Water_Amsol71_SM542R_Keywords)

    # write the AMSOL7.1 keywords for a SM5.42R point calculation in hexadecane to the AMSOL7.1 hexadecane input-file:
    Hexadecane_Amsol71_SM542R_Keywords = (
        """CHARGE=%s AM1 1SCF TLIMIT=15 GEO-OK SM5.42R\n& SOLVNT=GENORG IOFR=1.4345 ALPHA=0.00 BETA=0.00 GAMMA=38.93\n& DIELEC=2.06 FACARB=0.00 FEHALO=0.00 DEV\n"""
        % netcharge
    )
    Actual_Amsol71_InputFile_Hexadecane.write(Hexadecane_Amsol71_SM542R_Keywords)

    # write the name of the currently treated protonated state of the molecule into the AMSOL7.1 file
    # plus the number of atoms in the molecule
    logger.debug(
        f"len(ZmatMOPAC_Data) = number of atoms in molecule : {len(ZmatMOPAC_Data)}"
    )
    NumberOfAtomsInMolecule = len(ZmatMOPAC_Data)
    Molecule_Name_NrAtoms = "%s %d\n" % (MoleculeName, NumberOfAtomsInMolecule)
    Actual_Amsol71_InputFile_Water.write(Molecule_Name_NrAtoms)
    Actual_Amsol71_InputFile_Hexadecane.write(Molecule_Name_NrAtoms)

    # write a blank line after the keywords block and the line showing the name of the protonated state of the molecule to the AMSOL7.1 input-files
    blank_line = "\n"
    Actual_Amsol71_InputFile_Water.write(blank_line)
    Actual_Amsol71_InputFile_Hexadecane.write(blank_line)

    # write the lines of the MOPAC Z-matrix to the AMSOL7.1 input-files
    for line_keys in ZmatMOPAC_Data:
        Actual_Amsol71_InputFile_Water.write(ZmatMOPAC_Data[line_keys])
        Actual_Amsol71_InputFile_Hexadecane.write(ZmatMOPAC_Data[line_keys])

    Actual_Amsol71_InputFile_Water.close()
    Actual_Amsol71_InputFile_Hexadecane.close()

    logger.debug("just finished the function create_amsol71_inputfile(). ")

    return


def make_amsol71_input(path_and_file_ZmatMOPAC, MoleculeName):
    logger.debug("just entering main program in make_amsol71_input.py: ")
    ZmatMOPAC_data = {}
    ZmatMOPAC_data = read_ZmatMOPAC(path_and_file_ZmatMOPAC)
    logger.debug("calling create_amsol71_inputfile(): ")
    create_amsol71_inputfile(path_and_file_ZmatMOPAC, MoleculeName, ZmatMOPAC_data)
    return


if __name__ == "__main__":
    if len(sys.argv) != 3:  # if no input
        logger.error(" make_amsol71_input.py needs a ZmatMOPAC file as input.")
        logger.error(
            " The ZmatMOPAC file must have been generated by 'obabel ... -o mopin ...' (MOPAC Internals)"
        )
        logger.error(" Make sure that obabel is installed on your computer !!!")
        sys.exit()

    path_and_file_ZmatMOPAC = sys.argv[1]
    MoleculeName = sys.argv[2]
    make_amsol71_input(path_and_file_ZmatMOPAC, MoleculeName)
