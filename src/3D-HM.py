# Copyright (c) 2022, Sabine Reisser
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.


import logging
import os
import argparse
import subprocess
from config import *
import shutil
import inputgen
import hm_vector
import output
import surface

logging.basicConfig(format='%(levelname)s:%(name)s: %(message)s',
                    level=20)
logger = logging.getLogger(f"3D-HM")


if __name__ == "__main__":
    '''
    Calculate electrostatic potential (ESP) on solvent accessible surface of 
    molecule.
    If input is a pdb file, a pqr file with charges and radii is first
    using the given input force field.
    Finally, the ESP on the surface is evaluated to calculate a 3D hydrophobic
    moment vector.
    '''

    parser = argparse.ArgumentParser(description='Calculate 3D Hydrophobic '
                                                 'Moment Vector',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o', '--output', type=str,
                        default='output', help='output name')
    parser.add_argument('--die', type=float, default=78.54,
                        help='dielectric constant for APBS electrostatics '
                             'calculation')
    parser.add_argument('input', type=argparse.FileType('r'),
                        help='Input pdb/pqr/seq file'
                        )
    parser.add_argument('--ff', type=str,
                         choices=['AMBER', 'CHARMM', 'PARSE', 'custom'],
                         default='custom',
                         help='forcefield used by pdb2pqr'
                        )
    parser.add_argument('--neutraln',
                        action="store_true",
                        default=False,
                        help="Make the N-terminus of a protein neutral "
                             "(pdb2pqr option)"
                        )
    parser.add_argument("--neutralc",
                        action="store_true",
                        default=False,
                        help="Make the C-terminus of a protein neutral "
                             "(pdb2pqr option)"
    )
    parser.add_argument("--nocleanup",
                        dest='cleanup',
                        action="store_false",
                        default=True,
                        help="Keep all intermediate files of the calculation"
    )

    args = parser.parse_args()


    assert args.die > 0, 'Dielectric constant cannot be negative'

    file_name = args.input.name
    file_extension = file_name.split('.')[-1]
    output_name = args.output
    output_dir = f'{output_name}_OUT'

    logger.info(f'Creating output directory {output_dir} and copying '
                f'input file there.')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    shutil.copyfile(file_name, os.path.join(output_dir, file_name))

    os.chdir(output_dir)

    if file_extension == 'pqr':
        pqr_file = file_name
        inputgen.generate_apbs_input(pqr_file, 'apbs.in')
    else:
        if file_extension == 'pdb':
            pdb_file = file_name

        else:
            # expect amino acid sequence
            sequence = inputgen.read_sequence(file_name)
            logger.info(f'Read sequence: {sequence}')

            pdb_file = inputgen.create_pdb(sequence, output_name)

        # Use pdb2pqr to create pqr file
        pqr_file = inputgen.run_pdb2pqr(pdb_file, output_name, args)

    xyzr_file = inputgen.pqr2xyzr(pqr_file, output_name)

    # use NanoShaper to calculate surface

    ## create parameter file for NanoShaper
    with open(f'{DIR_3DHM}/templates/nanoshaper.prm', 'r') as template:
        temp = template.read()
        prm = temp.replace('XYZRFILE', xyzr_file)
    prm_file = f'{output_name}.prm'
    with open(prm_file, 'w') as prmfile:
        prmfile.write(prm)
    logger.info(f'Created input file {prm_file} for NanoShaper')

    ## run NanoShaper

    with open(f'{output_name}.NanoShaper', 'w') as out:
        return_code = subprocess.run([BIN_NANOSHAPER, prm_file],
                                  stdout=out,
                                  stderr=subprocess.STDOUT
                                  )


    # convert surface
    surface_file = 'triangulatedSurf.off'


    converted_surface_file = surface.write_converted_surface(surface_file, output_name)

    # set dielectric constant for APBS electrostatic potential calculation
    with open('apbs.in', 'r') as f:
        apbs_input = f.read()
        apbs_input = apbs_input.replace('78.5400', str(args.die))
    with open('apbs.in', 'w') as f:
        f.write(apbs_input)

    # run APBS
    with open(f'{output_name}.apbs', 'w') as out:
        return_code = subprocess.call([BIN_APBS, 'apbs.in'],
                                      stdout=out,
                                      stderr=subprocess.STDOUT
                                      )


    # project electrostatic potential onto surface
    multivalue_output = f'{output_name}.list_pot'
    return_code = subprocess.call([BIN_MULTIVALUE,
                                   converted_surface_file,
                                   f'{pqr_file}.dx',
                                   multivalue_output
                                   ],
                                  stdout=subprocess.DEVNULL)

    HM_vector, geometric_center = hm_vector.calc_hm_vector(multivalue_output,
                                                           output_name)

    output.write_pqr_with_hm(pqr_file,
                             HM_vector,
                             geometric_center,
                             output_name)

    if args.cleanup:
        logger.info('Removing intermediate files.')
        files_to_remove = [
            'apbs.in',
            'exposedIndices.txt',
            'exposed.xyz',
            'io.mc',
            'leap.log',
            f'{file_name}.dx',
            f'{output_name}.apbs',
            f'{output_name}.list_pot',
            f'{output_name}.NanoShaper',
            f'{output_name}.pqr.dx',
            f'{output_name}.prm',
            f'{output_name}.tleap_out',
            f'{output_name}.tri.list',
            f'{output_name}.xyzr',
            'stderror.txt',
            'tleap.in',
            'triangleAreas.txt',
            'triangulatedSurf.off',
        ]

        for f in files_to_remove:
            try:
                os.remove(f)
            except OSError:
                pass


