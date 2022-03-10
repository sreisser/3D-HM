import pdb2pqr
import re
import numpy as np
from config import *
import subprocess
import logging

logger = logging.getLogger(__name__)

def generate_3letter_code(seq, pH=7., nterm='N', cterm='C'):

    aa_map = {
        "A": "ALA",
        "C": "CYS",
        "D": "ASP",
        "E": "GLU",
        "F": "PHE",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "K": "LYS",
        "L": "LEU",
        "M": "MET",
        "N": "ASN",
        "P": "PRO",
        "Q": "GLN",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "U": "SEC",
        "V": "VAL",
        "W": "TRP",
        "Y": "TYR"
    }

    aa_ph_dependent = {
        # R
        "K" : [ 12.1, 'LYS', 'LYN'],
        "H" : [ 6.04, 'HIS', 'HID'],
        "E" : [ 4.15, 'GLH', 'GLU'],
        "D" : [ 3.71, 'ASH', 'ASP'],
        "C" : [ 8.14, 'CYS', 'CYM']
        # Y
    }

    three_letter_seq = []
    for letter in seq:
        if letter in aa_ph_dependent:
            if pH < aa_ph_dependent[letter][0]:
                three_letter_seq.append(aa_ph_dependent[letter][1])
            else:
                three_letter_seq.append(aa_ph_dependent[letter][2])
        else:
            three_letter_seq.append(aa_map[letter])

    seq_numbers = np.arange(1, len(seq)).astype(str)

    # take care of termini
    three_letter_seq[0] = f'{nterm}{three_letter_seq[0]}'
    three_letter_seq[-1] = f'{cterm}{three_letter_seq[-1]}'

    return " ".join(three_letter_seq), " ".join(seq_numbers)

def read_sequence(fname):

    with open(fname, 'r') as f:
        content = f.read()
        result = re.findall('^[A-Za-z\s\n]+$', content, flags=re.MULTILINE)

    sequence = ''.join(result[0].split()).upper()
    if not sequence:
        raise BaseException(f'Could not find valid sequence in {fname}')

    non_canonical_amino_acids = re.compile('.*[BJXO].*')
    if non_canonical_amino_acids.match(sequence):
        raise BaseException('Found non-canonical amino acid (BJXO)'
                            ' in sequence. Aborting.')
    return sequence

def create_pdb(sequence, output):
    three_letter_seq, seq_numbers = generate_3letter_code(sequence)
    output_pdb = f'{output}.pdb'


    with open(f'{DIR_3DHM}/templates/tleap_load.txt', 'r') as template:
        temp = template.read()
        tlp = temp.replace('SEQUENCE', three_letter_seq)
        tlp = tlp.replace('SEQ_NUMBERS', seq_numbers)
        tlp = tlp.replace('OUTPUT', output_pdb)
    tlp_file = 'tleap.in'
    with open(tlp_file, 'w') as f:
        f.write(tlp)

    with open(f'{output}.tleap_out', 'w') as out:
        return_code = subprocess.run([BIN_TLEAP, '-f', tlp_file],
                                      stdout=out,
                                      stderr=subprocess.STDOUT
                                      )
        return output_pdb

def generate_apbs_input(pqr, output):
    size_obj = pdb2pqr.psize.Psize()
    size_obj.run_psize(pqr)
    input = pdb2pqr.inputgen.Input(pqr, size_obj, 'mg-auto', False, potdx=True)
    input.print_input_files(output)



def run_pdb2pqr(pdb_file, output, args):
    pqr_file = f'{output}.pqr'

    pdb2pqr_parser = pdb2pqr.main.build_main_parser()
    pdb2pqr_args = ['--apbs-input=apbs.in']

    if not args.ff or args.ff == 'custom':
        pdb2pqr_args += [
                f'--userff={DIR_3DHM}/dat/{args.ff}.DAT',
                f'--usernames={DIR_3DHM}/dat/{args.ff}.names'
        ]
    else:
        pdb2pqr_args += [f'--ff={args.ff}']

    if args.neutraln:
        if protonated(pdb_file):
            logger.warning('File is already protonated, cannot change N-terminus')
        else:
            pdb2pqr_args.append('--neutraln')
    if args.neutralc:
        if protonated(pdb_file):
            logger.warning('File is already protonated, '
                         'cannot change C-terminus')
        else:
            pdb2pqr_args.append('--neutralc')

    params = pdb2pqr_parser.parse_args(
        pdb2pqr_args +
        [
            pdb_file,
            pqr_file
        ]
    )

    print(pdb2pqr_args)

    # Loading topology files
    definition = pdb2pqr.io.get_definitions()

    pdblist, is_cif = pdb2pqr.io.get_molecule(pdb_file)
    # drop water
    pdblist = pdb2pqr.main.drop_water(pdblist)
    # Setting up molecule
    biomolecule, definition, ligand = pdb2pqr.main.setup_molecule(
        pdblist, definition, None
    )

    # Setting termini states for biomolecule chains
    biomolecule.set_termini(params.neutraln, params.neutralc)

    results = pdb2pqr.main.non_trivial(
        args=params,
        biomolecule=biomolecule,
        ligand=ligand,
        definition=definition,
        is_cif=is_cif,
    )

    pdb2pqr.main.print_pqr(
        args=params,
        pqr_lines=results["lines"],
        header_lines=results["header"],
        missing_lines=results["missed_residues"],
        is_cif=is_cif,
    )

    if params.apbs_input:
        pdb2pqr.io.dump_apbs(params.output_pqr, params.apbs_input)

  #  logging.basicConfig(level=20)
    return pqr_file


def pqr2xyzr(pqr_file, output):

    pqrfile_handle = open(pqr_file, "r")
    pqrfile_content = pqrfile_handle.readlines()
    pqrfile_handle.close()

    xyzr_name = f'{output}.xyzr'

    re_pqr = re.compile(
        "^ATOM\s{2}([0-9\s]{5})\s([A-Z0-9\s]{4}).([A-Z\s]{4}).([\-0-9\s]{4})"
        ".\s{3}([0-9\-\.\s]{8})([0-9\-\.\s]{8})([0-9\-\.\s]{8})\s+([0-9\.\-]+)"
        "\s+([0-9\.\-]+)\s*$")

    total_charge = 0.
    xyzrfile_content = ""

    res_charge = 0
    old_res = 0
    for line in pqrfile_content:
        atom = re_pqr.match(line)
        if atom:
            x = float(atom.group(5))
            y = float(atom.group(6))
            z = float(atom.group(7))
            radius = float(atom.group(9))
            xyzrfile_content += "%f %f %f %f\n" % (x, y, z, radius)

            # check charge
            total_charge += float(atom.group(8))

            res = atom.group(4)
            if res != old_res:
                print(f'{atom.group(4)}, {np.round(res_charge, 4)}')
                res_charge = float(atom.group(8))
            else:
                res_charge += float(atom.group(8))
            old_res = res


    if np.abs(np.round(total_charge, 6)) % 1 != 0:
        logger.warning(f'Total charge is not integer: {total_charge}!')
    else:
        logger.info(f'Total charge: {total_charge:.2f}')

    xyzr_file = open(xyzr_name, "w")
    xyzr_file.write(xyzrfile_content)
    xyzr_file.close()
    logger.info(f"xyzr file generated from file {pqr_file}")

    return xyzr_name


def protonated(pdb_file):
    with open(pdb_file, 'r') as f:
        content = f.read()
    re_hydrogen = re.compile(
        "ATOM\s{2}([0-9\s]{5})\s(H[A-Z0-9\s]{3}|\sH[A-Z0-9\s]{2}).*")
    result = re_hydrogen.findall(content, re.MULTILINE)

    if result:
        return True
    else:
        return False

