#!/usr/bin/env python
import argparse
import os
import shutil
import sys
import logging
logger = logging.getLogger("prepare")  # Placeholder for future betterness

def prepare(molfile, namedata=None, namepath=None, tempdir=os.getcwd()):
    tempdir = tempdir.rstrip("/")

    # Create working directory
    if not os.path.isdir(tempdir):
        logger.debug("Creating molecule temp dir: {}".format(tempdir))
        os.makedirs(tempdir)

    srcbase = os.path.basename(molfile)
    tempmol = os.path.join(tempdir, srcbase)
    namefile = os.path.join(tempdir, "name.txt")

    # Bring mol2 into working directory
    if not os.path.exists(tempmol) or not os.path.samefile(molfile, tempmol):
        logger.debug("Copying mol file to temp dir: {}".format(tempmol))
        shutil.copy(molfile, tempmol)

    if namedata is not None:
        logger.debug(">>>>>> Writing name.txt")
        f = open(namefile, "w")
        f.write(namedata + "\n")
        f.close()
    elif namepath is not None and not os.path.samefile(namepath, namefile):
        logger.debug("Copying namefile from {0}".format(namepath))
        shutil.copy(namepath, namefile)
    else:
        logger.warning("Skipping name.txt generation")


def create_namedata(name, smiles, longname, prot_id=None, temp_id=None):
    if temp_id is None:
        temp_id = 1
    if prot_id is None:
        line = f"name.txt {temp_id} {name} {smiles} | {longname}"
    else:
        line = f"name.cxcalc.txt {temp_id} {name} {prot_id} {smiles} _ | {longname}"
    return line


def main(args):
    parser = argparse.ArgumentParser("Prepare mol2 for database generation")
    parser.add_argument("mol2", help="Path to mol2 file")
    parser.add_argument(
        "-d",
        "--dir",
        nargs="?",
        default=os.getcwd(),
        help="Directory to prpare ligand in",
    )
    parser.add_argument(
        "-f", "--namefile", default=None, help="Existing namefile to use"
    )
    parser.add_argument(
        "-S",
        "--smilesfile",
        default=None,
        help="Read name and SMILES from SMILES file (first line)",
    )
    parser.add_argument("-n", "--name", default=None, help="Ligand name")
    parser.add_argument("-p", "--prot_id", default=None, help="Ligand protomer ID")
    parser.add_argument("-s", "--smiles", default=None, help="Ligand SMILES")
    parser.add_argument(
        "-l", "--long_name", default="NO_LONG_NAME", help="Ligand long name"
    )
    parser.add_argument(
        "-t", "--temp_id", default=None, help="Ligand temporary (internal) ID"
    )
    params = parser.parse_args(args)

    if params.namefile is None:
        if params.smilesfile is not None:
            sf = open(params.smilesfile)
            smiles, name = next(sf).split()[0:2]
        elif params.name is None:
            logger.error("No name or namefile provided!")
            sys.exit(-1)
        elif params.smiles is None:
            logger.warn("SHOULD include ligand smiles OR namefile. Using empty smiles")
            smiles = "none"
        else:
            name = params.name
            smiles = params.smiles
        namefile = None
        namedata = create_namedata(
            name=name,
            smiles=smiles,
            longname=params.long_name,
            prot_id=params.prot_id,
            temp_id=params.temp_id,
        )
    else:
        namefile = params.namefile
        namedata = None
    prepare(params.mol2, namedata=namedata, namepath=namefile, tempdir=params.dir)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
