"""
This is a python script that converts u(rho, T), P(rho, T), Cs(rho,T), S(rho, T)
to T(rho, u), P(rho, u), Cs(rho, u), S(rho, u), which is more useful for SPH calculations
"""


import matplotlib.pyplot as plt
from collections import OrderedDict
import numpy as np
import pandas as pd
import csv
import sys
from scipy.interpolate import interp1d
from scipy import interpolate


def emptyLineIndices(f):
    empty_lines = [0]
    with open(f, 'r') as infile:
        reader = csv.reader(infile)
        next(reader)  # drop header row
        for index, row in enumerate(reader):
            if len(row) == 0:
                empty_lines.append(index)
    infile.close()
    return empty_lines


def chunkFile(f, emtpy_lines):
    densities = []
    d = {}
    with open(f, 'r') as infile:
        reader = csv.reader(infile)
        headers = next(reader)
        reader = list(reader)
        for index, line in enumerate(empty_lines):

            temp_dict = {}
            for i in headers:
                temp_dict.update({i: []})

            if (index + 1) != len(empty_lines):
                min, max = empty_lines[index] + 1, empty_lines[index + 1] - 1
                trimmed_reader = reader[min:max]
                for row in trimmed_reader:
                    for index2, i in enumerate(row):
                        header = headers[index2]
                        temp_dict[header].append(reformat(i))
                density = reformat(temp_dict['Pressure (Pa)'][0])
                densities.append(density)
                d.update({density: temp_dict})




    return d


def reformat(number):
    if isinstance(number, str):
        if '-101' in str(number):
            new_num  = float(number.split('-')[0]) * (10**(-101))
            return new_num
        else:
            return float(number)
    else:
        return number


def recalculateEnergies(d, grid_number, min_energy, delta):
    """
    For each density sample, we want the same exponential energy grid
    :param d:
    :param grid_number:
    :param min_energy:
    :param delta:
    :return:
    """
    densities = d.keys()
    new_energies = []
    for i in range(0, grid_number):
        new_energy = min_energy * (delta**i)
        new_energies.append(new_energy)
    for i in densities:
        d[i].update({'Energy (J/kg)': new_energies})

    return d




nu = 120  # number of the grid for the internal energy (exponential)

infile_path = 'granite.table.csv'
empty_lines = emptyLineIndices(f=infile_path)
sorted_dict = chunkFile(f=infile_path, emtpy_lines=empty_lines)
densities = sorted_dict.keys()

infile_df = pd.read_csv(infile_path)
energy = [reformat(i) for i in list(infile_df['Energy (J/kg)'])]
min_energy = min(energy)
max_energy = max(energy)
delta = (min_energy / max_energy)**(1/(nu-1))

sorted_dict = recalculateEnergies(d=sorted_dict, grid_number=nu, min_energy=min_energy, delta=delta)

for i in densities:

    energies = sorted_dict[i]['Energy (J/kg)']
    temperatures = sorted_dict[i]['Temperature (K)']
    pressures = sorted_dict[i]['Pressure (Pa)']
    sound_speeds = sorted_dict[i]['Sound speed (m/s)']
    entropies = sorted_dict[i]['Entropy (J/kg/K)']


    f_temperature = interpolate.interp1d(energies, temperatures, kind='linear', fill_value='extrapolate')
    sorted_dict[i].update({'Temperature (K)': f_temperature(energies)})


    f_pressure = interpolate.interp1d(temperatures, pressures, kind='linear', fill_value='extrapolate')
    sorted_dict[i].update({'Pressure (Pa)': f_pressure(sorted_dict[i]['Temperature (K)'])})


    f_soundspeed = interpolate.interp1d(temperatures, sound_speeds, kind='linear', fill_value='extrapolate')
    sorted_dict[i].update({'Sound speed (m/s)': f_soundspeed(sorted_dict[i]['Temperature (K)'])})

    f_entropy = interpolate.interp1d(temperatures, entropies, kind='linear', fill_value='extrapolate')
    sorted_dict[i].update({'Entropy (J/kg/K)': f_entropy(sorted_dict[i]['Temperature (K)'])})








# infile_df = pd.read_csv(infile_path)
#
# density = sorted(list(set([reformat(i) for i in list(infile_df['Density (kg/m3)'])])))  # remove duplicates, then sort
# temperature = sorted(list(set([reformat(i) for i in list(infile_df['Temperature (K)'])])))
# energy = [reformat(i) for i in list(infile_df['Energy (J/kg)'])]
# pressure = [reformat(i) for i in list(infile_df['Pressure (Pa)'])]
# sound_speed = [reformat(i) for i in list(infile_df['Sound speed (m/s)'])]
# entropy = [reformat(i) for i in list(infile_df['Entropy (J/kg/K)'])]
#
# min_energy = min(energy)
# max_energy = max(energy)
# delta = (min_energy / max_energy)**(1 / (nu - 1))
#
# new_energy = [min_energy * (delta**i) for i in range(0, nu)]
#
# new_temperature = []
# new_pressure = []
# new_sound_speed = []
# new_entropy = []
#
# for m in range(0, nu):
#
#     # internal energy
#     f_temperature = interpolate.interp1d(energy[m:], temperature[m:], kind='linear', fill_value='extrapolate')
#     new_temperature.append(f_temperature(new_energy))
#
#     # pressure
#     f_pressure = interpolate.interp1d(temperature[m:], pressure[m:], kind='linear', fill_value='extrapolate')
#     new_pressure.append(f_pressure(new_temperature[m]))
#
#     # sound speed
#     f_soundspeed = interpolate.interp1d(temperature[m:], sound_speed[m:], kind='linear', fill_value='extrapolate')
#     new_sound_speed.append(f_soundspeed(new_temperature[m]))
#
#     # entropy
#     f_entropy = interpolate.interp1d(temperature[m:], entropy[m:], kind='linear', fill_value='extrapolate')
#     new_entropy.append(f_entropy(new_temperature[m]))
#
# new_temperature = np.array(new_temperature)
# new_pressure = np.array(new_pressure)
# new_sound_speed = np.array(new_sound_speed)
# new_entropy = np.array(new_entropy)
#
# for m in range(0, len(density), int(len(density)/6)):
#
#     ax = [0, 0, 0, 0]
#
#     fig = plt.figure(figsize = (10,6.128))
#
#     ax[0] = fig.add_subplot(221)
#     ax[1] = fig.add_subplot(222)
#     ax[2] = fig.add_subplot(223)
#     ax[3] = fig.add_subplot(224)
#
#     ax[0].semilogy(np.array(temperature) * 1e-3, np.array(energy[m:]) * 1e-6, '--', label="original ANEOS")
#     ax[0].semilogy(new_temperature[m:] * 1e-3, np.array(new_energy[m:]) * 1e-6, '-.', label="modified")
#     ax[1].semilogy(np.array(temperature) * 1e-3, np.array(pressure[m:]) * 1e-6,'--', new_temperature[m:] * 1e-3, new_pressure[m:] * 1e-6,'-.')
#     ax[2].plot(np.array(temperature) * 1e-3, np.array(sound_speed[m:]) * 1e-3,'--', new_temperature[m:] * 1e-3, new_sound_speed[m:] * 1e-3,'-.')
#     ax[3].plot(np.array(temperature) * 1e-3, np.array(entropy[m:]) * 1e-3,'--', new_temperature[m:] * 1e-3, new_entropy[m:] * 1e-3,'-.')
#
#     ax[0].legend(frameon=False)
#
#     ax[0].set_ylabel('Energy (MJ/kg)', fontsize=10)
#     ax[1].set_ylabel('Pressure (MPa)', fontsize=10)
#     ax[2].set_ylabel('Sound Speed (km/s)', fontsize=10)
#     ax[3].set_ylabel('Entropy (kJ/K/kg)', fontsize=10)
#     ax[2].set_xlabel('Temperature ($10^3$ K)', fontsize=10)
#     ax[3].set_xlabel('Temperature ($10^3$ K)',fontsize=10)
#
#     fig.suptitle("Density: %3.3f kg/m$^3$" %(density[m]))
#     # plt.show()
#     # fig.savefig("Density" + str(m) + ".png")
