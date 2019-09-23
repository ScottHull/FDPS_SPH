import os
import pandas as pd
from math import sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

interpolation_file = pd.read_csv('granite.rho_u.csv')


def bilinear_interpolation(x, y, points):

    '''
    Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0
    '''

    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
            ) / ((x2 - x1) * (y2 - y1) + 0.0)


def weights(x, y, points):

    # R1 = points[0][2] + (((x - points[0][0]) / (points[1][0] - points[0][0])) * (points[1][2] - points[0][2]))
    # R2 = points[2][2] + (((x - points[0][0]) / (points[1][0] - points[0][0])) * (points[3][2] - points[2][2]))
    print(points)
    R1 = (((points[1][0] - x) / (points[1][0] - points[0][0])) * points[0][2]) + (((x - points[0][0]) / (points[1][0] - points[0][0])) * points[1][2])
    R2 = (((points[1][0] - x) / (points[1][0] - points[0][0])) * points[2][2]) + (((x - points[0][0]) / (points[1][0] - points[0][0])) * points[3][2])

    return R1, R2



def matrix_variable(density_array, energy_array, variable_array):


    m = []
    row = []
    count = 0
    density_index = 0
    energy_index = 0
    current_density = None
    for index, i in enumerate(density_array):
        if count == 0:
            current_density = i
        if i != current_density:
            current_density = i
            density_index += 1
            energy_index = 0
            m.append(row)
            row = []

        row.append(variable_array[count])

        energy_index += 1
        count += 1

    return m


def get_database_boundaries(density_array, energy_array, variable_matrix):

    density_neighbors = (density_array[0], density_array[len(density_array) - 1])
    energy_neighbors = (energy_array[0], energy_array[len(energy_array) - 1])

    corresponding_variables = (variable_matrix[0][0], variable_matrix[len(variable_matrix) - 1][0],
                                variable_matrix[0][60 - 1],
                                variable_matrix[len(variable_matrix) - 1][60 - 1])

    q11 = (density_neighbors[0], energy_neighbors[0], corresponding_variables[0])
    q12 = (density_neighbors[1], energy_neighbors[0], corresponding_variables[1])
    q21 = (density_neighbors[0], energy_neighbors[1], corresponding_variables[2])
    q22 = (density_neighbors[1], energy_neighbors[1], corresponding_variables[3])
    return [q11, q12, q21, q22]


def plot_interpolate(points, interpolated_point, weights):

    fig = plt.figure()
    ax = Axes3D(fig)


def plot_database(densities, internal_energies, variables):

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot(densities, internal_energies, variables)



def bilinear_interpolate(density, energy, density_array, energy_array, variable_array):

    interpolated_point = (density, energy)

    variable_matrix = matrix_variable(density_array=density_array, energy_array=energy_array,
                                        variable_array=variable_array)


    points = get_database_boundaries(
        density_array=density_array, energy_array=energy_array, variable_matrix=variable_matrix)

    w = weights(density, energy, points)
    i = bilinear_interpolation(density, energy, points)


    plot_database(densities=density_array, internal_energies=energy_array, variables=variable_array)

    plt.show()





density_array = interpolation_file['Density (kg/m3)']
energy_array = interpolation_file['Energy (J/kg)']
temperature_array = interpolation_file['Temperature (K)']
bilinear_interpolate(density=2000, energy=1.25e7, density_array=density_array, energy_array=energy_array,
                     variable_array=temperature_array)