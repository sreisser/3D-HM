# Copyright (c) 2022, Sabine Reisser
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.


import numpy as np
import re


def write_pqr_with_hm(pqr_file, hm_vector, geometric_center, output_name):

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
