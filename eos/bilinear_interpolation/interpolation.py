import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class BilinearInterpolation:

    def __init__(self, density_array, internal_energy_array, variable_array, density, internal_energy):

        self.density_array = density_array
        self.internal_energy_array = internal_energy_array
        self.variable_array = list(variable_array)
        self.density = density
        self.internal_energy = internal_energy

        self.variable_matrix = self.matrix_variable()
        self.points = self.getBoundaries()


    def matrix_variable(self):

        m = []
        row = []
        count = 0
        density_index = 0
        energy_index = 0
        current_density = None
        for index, i in enumerate(self.density_array):
            if count == 0:
                current_density = i
            if i != current_density:
                current_density = i
                density_index += 1
                energy_index = 0
                m.append(row)
                row = []

            row.append(self.variable_array[count])

            energy_index += 1
            count += 1

        return m

    def getBoundaries(self):
        
        density_neighbors = (self.density_array[0], self.density_array[len(self.density_array) - 1])
        energy_neighbors = (self.internal_energy_array[0], self.internal_energy_array[len(self.internal_energy_array) - 1])
    
        corresponding_variables = (self.variable_matrix[0][0], self.variable_matrix[len(self.variable_matrix) - 1][0],
                                   self.variable_matrix[0][60 - 1],
                                   self.variable_matrix[len(self.variable_matrix) - 1][60 - 1])
    
        q11 = (density_neighbors[0], energy_neighbors[0], corresponding_variables[0])
        q12 = (density_neighbors[1], energy_neighbors[0], corresponding_variables[1])
        q21 = (density_neighbors[0], energy_neighbors[1], corresponding_variables[2])
        q22 = (density_neighbors[1], energy_neighbors[1], corresponding_variables[3])
        return [q11, q12, q21, q22]

    def interpolate(self):

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

        points = sorted(self.points)               # order points by x, then by y
        (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

        if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
            raise ValueError('points do not form a rectangle')
        if not x1 <= self.density <= x2 or not y1 <= self.internal_energy <= y2:
            raise ValueError('(x, y) not within the rectangle')

        return (q11 * (x2 - self.density) * (y2 - self.internal_energy) +
                q21 * (self.density - x1) * (y2 - self.internal_energy) +
                q12 * (x2 - self.density) * (self.internal_energy - y1) +
                q22 * (self.density - x1) * (self.internal_energy - y1)
                ) / ((x2 - x1) * (y2 - y1) + 0.0)







interpolation_file = pd.read_csv('/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/FDPS_SPH/eos/granite.rho_u.csv')
original_data_file = pd.read_csv('/Users/scotthull/Documents - Scott’s MacBook Pro/PhD Research/FDPS_SPH/eos/granite.table.csv')

density_array = interpolation_file['Density (kg/m3)']
energy_array = interpolation_file['Energy (J/kg)']
temperature_array = interpolation_file['Temperature (K)']

# density = 1520
# internal_energy = 1e9

# b = BilinearInterpolation(density=density, internal_energy=internal_energy, density_array=density_array,
#                           internal_energy_array=energy_array, variable_array=temperature_array)

fig = plt.figure()
ax = Axes3D(fig)
ax.plot(density_array, energy_array, temperature_array, color='blue')
ax.plot(original_data_file['Density (kg/m3)'], original_data_file['Energy (J/kg)'],
           original_data_file['Temperature (K)'], color='red')
for index, i in enumerate(density_array):
    density = density_array[index]
    internal_energy = energy_array[index]
    b = BilinearInterpolation(density=density, internal_energy=internal_energy, density_array=density_array,
                              internal_energy_array=energy_array, variable_array=temperature_array)
    ax.scatter(density, internal_energy, b.interpolate(), color='green')
ax.set_xlabel("Density (kg/m3)")
ax.set_ylabel("Internal Energy (J/kg)")
ax.set_zlabel("Temperature (K)")

plt.show()





