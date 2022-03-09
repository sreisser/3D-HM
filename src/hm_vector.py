import numpy as np
import re
import datetime
import logging

logger = logging.getLogger(__name__)

def triangle_area(triangle):
    a = np.linalg.norm(triangle[1, :3] - triangle[0, :3])
    b = np.linalg.norm(triangle[2, :3] - triangle[0, :3])
    c = np.linalg.norm(triangle[2, :3] - triangle[1, :3])
    s = 0.5 * (a + b + c)
    return np.sqrt(s * (s - a) * (s - b) * (s - c))


def calc_hm_vector(multivalue_output, output_name):
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
    output += f'\tAverage absolute ESP on surface: {average_esp_abs:.3f}\n'
    output += f'\tGeometric molecule center: ' \
              f'{" ".join(geometric_center_str)}\n'
    output += f'\tHydrophobic moment vector: ' \
              f'{" ".join(hm_vector_str)}\n'
    output += f'\tAbsolute HM vector {abs_hm_vector:.3f} AkT/e\n'
    output += f'\tSurface: {n_triangles} triangles, ' \
              f'area: {total_surface:.3f} A**2\n'
    output += f'\tAngle between HM vector and z-axis: {angle_z_axis:.3f} degrees\n'

    with open(f'HM_{output_name}.output', 'w') as f:
        f.write(output)

    logger.info(output)

    return hydrophobic_moment_vector, geometric_center


def angle_with_z(hm_vector):
    z_axis = np.array([0, 0, 1])
    return np.arccos(np.dot(hm_vector, z_axis) /
                     (np.linalg.norm(hm_vector) * np.linalg.norm(z_axis))) \
           * 180 / np.pi