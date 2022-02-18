# import requests
import pdb2pqr
import argparse
import re
import subprocess
import numpy as np
import datetime

def check_filename(fname):
    choices = ['pdb', 'pqr']
    ext = fname.name.split('.')[-1]
    if ext not in choices:
        return False
    return ext


def pqr2xyzr(pqr_file):
    pqrfile_handle = open(pqr_file, "r")
    pqrfile_content = pqrfile_handle.readlines()
    pqrfile_handle.close()

    xyzr_name = pqr_file.replace('.pqr', '.xyzr')

    re_pqr = re.compile(
        "^ATOM\s{2}([0-9\s]{5})\s([A-Z0-9\s]{4}).([A-Z\s]{4}).([0-9\s]{4})"
        ".\s{3}([0-9\-\.\s]{8})([0-9\-\.\s]{8})([0-9\-\.\s]{8})\s+([0-9\.\-]+)"
        "\s+([0-9\.\-]+)\s*$")

    xyzrfile_content = ""
    for line in pqrfile_content:
        atom = re_pqr.match(line)
        if atom:
            x = float(atom.group(5))
            y = float(atom.group(6))
            z = float(atom.group(7))
            radius = float(atom.group(9))

            xyzrfile_content += "%f %f %f %f\n" % (x, y, z, radius)

    xyzr_file = open(xyzr_name, "w")
    xyzr_file.write(xyzrfile_content)
    xyzr_file.close()
    print(f"xyzr file generated from file {pqr_file}")

    return xyzr_name


def run_pdb2pqr(pdb_file, pqr_file):
    pdb2pqr_parser = pdb2pqr.main.build_main_parser()

    if not args.ff or args.ff == 'custom':
        args_ff = [
                f'--userff=../dat/{args.ff}.DAT',
                f'--usernames=../dat/{args.ff}.names'
        ]
    else:
        args_ff = [f'--ff={args.ff}']

    params = pdb2pqr_parser.parse_args(
        [
            '--assign-only',
            '--apbs-input=apbs.in'
        ] +
        args_ff +
        [
            pdb_file,
            pqr_file
        ]
    )
    print(params)

    # Loading topology files
    definition = pdb2pqr.io.get_definitions()

    pdblist, is_cif = pdb2pqr.io.get_molecule(file_name)
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

    #  print(results['header'])

    pdb2pqr.main.print_pqr(
        args=params,
        pqr_lines=results["lines"],
        header_lines=results["header"],
        missing_lines=results["missed_residues"],
        is_cif=is_cif,
    )

    if params.apbs_input:
        pdb2pqr.io.dump_apbs(params.output_pqr, params.apbs_input)


def write_converted_surface(surface_file, converted_surface_file):
    with open(surface_file, 'r') as f:
        surface = f.read()

    reg_vert = re.compile(
        "(\-*[0-9]+\.[0-9]+)\s+(\-*[0-9]+\.[0-9]+)\s+(\-*[0-9]+\.[0-9]+)\s*\n")
    reg_tri = re.compile("3\s([0-9]+)\s+([0-9]+)\s+([0-9]+)\s*\n")

    vertices = reg_vert.findall(surface)
    vertices = np.array(vertices)

    triangles = reg_tri.findall(surface)
    triangles = np.array(triangles).astype(int)

    output = ''
    for triangle in triangles:
        for i in triangle:
            output += ",".join(vertices[i]) + '\n'

    with open(converted_surface_file, 'w') as f:
        f.write(output)

def calc_hm_vector(multivalue_output):
    with open(multivalue_output, 'r') as f:
        esp_on_surface = f.read()

    reg_coordinates = re.compile(
        "([0-9\.\-e\+]+),([0-9\.\-e\+]+),([0-9\.\-e\+]+),([0-9\.\-e\+]+)\n")

    coordinates_esp = reg_coordinates.findall(esp_on_surface)
    # reshaping generates array of triangles:
    # 3 corners form one triangle
    # one corner has 4 values: x, y, z, ESP
    coordinates_esp = np.array(coordinates_esp).astype(float) \
        .reshape((-1, 3, 4))

    # need these three arrays for HM vector calculation
    esp_absolute = np.abs(coordinates_esp[:, :, 3]).reshape(-1, 3, 1)
    coordinates = coordinates_esp[:, :, :3]
    triangle_areas = np.array(list(map(triangle_area, coordinates)))
    triangle_areas = triangle_areas.reshape(-1, 1, 1)
    n_triangles = triangle_areas.shape[0]

    total_surface = np.sum(triangle_areas)
    sum_vertex_vectors = np.sum(coordinates * triangle_areas / 3, axis=(1, 0))
    geometric_center = sum_vertex_vectors / total_surface
    average_esp_abs = np.sum(esp_absolute * triangle_areas / 3) \
                      / total_surface

    sum_hm_vectors = np.sum(
        coordinates * esp_absolute * triangle_areas / 3, axis=(1, 0))

    hydrophobic_moment_vector = -(sum_hm_vectors -
                                  average_esp_abs * sum_vertex_vectors) / \
                                total_surface
    abs_hm_vector = np.linalg.norm(hydrophobic_moment_vector)

    angle_z_axis = angle_with_z(hydrophobic_moment_vector)
    geometric_center_str = [f'{i:.3f}' for i in geometric_center]
    hm_vector_str = [f'{i:.3f}' for i in hydrophobic_moment_vector]

    output = f'3D hydrophobic moment vector calculation performed on ' \
             f'{datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")}\n'
    output += f'Average absolute ESP on surface: {average_esp_abs:.3f}\n'
    output += f'Geometric molecule center: ' \
              f'{" ".join(geometric_center_str)}\n'
    output += f'Hydrophobic moment vector: ' \
              f'{" ".join(hm_vector_str)}\n'
    output += f'Absolute HM vector {abs_hm_vector:.3f} AkT/e\n'
    output += f'Surface: {n_triangles} triangles, ' \
              f'area: {total_surface:.3f} A**2\n'
    output += f'Angle between HM vector and z-axis: {angle_z_axis:.3f} degrees'

    with open(f'HM_{output_name}.output', 'w') as f:
        f.write(output)

    print(output)



    return hydrophobic_moment_vector, geometric_center


def triangle_area(triangle):
    a = np.linalg.norm(triangle[1, :3] - triangle[0, :3])
    b = np.linalg.norm(triangle[2, :3] - triangle[0, :3])
    c = np.linalg.norm(triangle[2, :3] - triangle[1, :3])
    s = 0.5 * (a + b + c)
    return np.sqrt(s * (s - a) * (s - b) * (s - c))

def angle_with_z(hm_vector):
    z_axis = np.array([0, 0, 1])
    return np.arccos(np.dot(hm_vector, z_axis) /
                   (np.linalg.norm(hm_vector) * np.linalg.norm(z_axis))) \
           * 180 / np.pi

def write_pqr_with_hm(pqr_file, hm_vector, geometric_center):

    with open(pqr_file, 'r') as f:
        pqr_content = f.read()

    # print all atoms
    reg_atom = re.compile("ATOM.*\n")
    atoms = reg_atom.findall(pqr_content)
    output = "".join(atoms)

    # print pseudo-atoms forming the HM vector
    hm_abs = np.linalg.norm(hm_vector)
    hm_normalized = hm_vector/hm_abs
    n_points = int(np.ceil(hm_abs)) + 1
    dist_points = hm_abs / (n_points - 1)

    for i in range(n_points):
        if i == 0:
            # vector starts in geometric center of molecule
            output += "ATOM    999  GMC HYM    99    " \
                      f"{geometric_center[0]:8.3f}" \
                      f"{geometric_center[1]:8.3f}" \
                      f"{geometric_center[2]:8.3f}  0.0000 1.0000\n"
        else:
            if i < n_points - 1:
                name = 'X'
            else:
                name = 'TIP'
            point = geometric_center + i * dist_points * hm_normalized
            output += f"ATOM    999  {name:3s} HYM    99    " \
                      f"{point[0]:8.3f}" \
                      f"{point[1]:8.3f}" \
                      f"{point[2]:8.3f}  0.0000 1.0000\n"

    with open(f'{output_name}_with_HM.pqr', 'w') as f:
        f.write(output)

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
    parser.add_argument('-o', '--output', metavar='output', type=str,
                        default='output', help='output name')
    parser.add_argument('--die', type=float, default=78.54,
                        help='dielectric constant')
    parser.add_argument('input', type=argparse.FileType('r'),
                        help='Input pdb file/pqr file')
    parser.add_argument('--ff', type=str,
                         choices=['AMBER', 'CHARMM', 'PARSE', 'custom'],
                         default='custom',
                         help='forcefield')

    args = parser.parse_args()
    print(args)

    assert args.die > 0, 'Dielectric constant cannot be negative'
    assert check_filename(args.input), 'Invalid filename (expects pqr or pdb)'

    file_name = args.input.name
    file_type = check_filename(args.input)
    output_name = args.output

    if file_type == 'pqr':
        pqr_file = file_name
    else:
        # Use pdb2pqr to create pqr file
        pqr_file = f'{output_name}.pqr'
        run_pdb2pqr(file_name, pqr_file)

    xyzr_file = pqr2xyzr(pqr_file)

    # use NanoShaper to calculate surface

    ## create parameter file for NanoShaper
    # TODO: get NanoShaper 0.7 config file
    with open('../doc/TEMPLATE.prm', 'r') as template:
        temp = template.read()
        prm = temp.replace('XYZRFILE', xyzr_file)
    prm_file = f'{output_name}.prm'
    with open(prm_file, 'w') as prmfile:
        prmfile.write(prm)
    print(f'Created input file {prm_file} for NanoShaper')

    ## run NanoShaper

    return_code = subprocess.call(['NanoShaper', prm_file])
    print("Output of call() : ", return_code)

    # convert surface
    surface_file = 'triangulatedSurf.off'

    converted_surface_file = f'{output_name}.tri.list'
    write_converted_surface(surface_file, converted_surface_file)

    # set dielectric constant for APBS electrostatic potential calculation
    with open('apbs.in', 'r') as f:
        apbs_input = f.read()
        apbs_input = apbs_input.replace('78.5400', str(args.die))
    with open('apbs.in', 'w') as f:
        f.write(apbs_input)

    # run APBS
    # TODO: set apbs_bin in config
    apbs_dir = '/home/sabine/3D-HM/APBS-3.4.0.Linux'
    apbs_bin = f'{apbs_dir}/bin/apbs'
    return_code = subprocess.call([apbs_bin, 'apbs.in'])
    print("Output of call() : ", return_code)

    multivalue_bin = f'{apbs_dir}/share/apbs/tools/bin/multivalue'
    multivalue_output = f'{output_name}.list_pot'
    return_code = subprocess.call([multivalue_bin,
                                   converted_surface_file,
                                   f'{pqr_file}.dx',
                                   multivalue_output
                                   ],
                                  stdout=subprocess.DEVNULL)

    hm_vector, geometric_center = calc_hm_vector(multivalue_output)

    write_pqr_with_hm(pqr_file, hm_vector, geometric_center)
