import numpy as np
import re

def write_converted_surface(surface_file, output_name):

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

    converted_surface_file = f'{output_name}.tri.list'
    with open(converted_surface_file, 'w') as f:
        f.write(output)

    return converted_surface_file



