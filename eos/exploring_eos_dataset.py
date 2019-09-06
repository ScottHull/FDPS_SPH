import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
from scipy.interpolate import interp1d
from scipy import interpolate

# ----- A user has to change these three parameters  ----------------

inputfilename = "granite.table.txt"  # input ANEOS file. This follows the format from iSALE
outputfilename = "granite.rho_u.txt"  # output ANEOS file
nu = 120  # number of the grid for the internal energy (exponential)

# This function is to correct the original ANEOS format that does not include "E"
# This seems to  occur when the exponent reaches -101
def reformat(number):
    """
    Fixes a bug in an input file.
    :param number:
    :return:
    """
    if number.find('E') == -1:
        exponent = "-101"
        mantissa = number.split(exponent)
        return float(mantissa[0])*10**float(exponent)
    else:
        mantissa, exponent =  number.split('E')

    return float(mantissa)*10**float(exponent)


aneosfile =  [line.split() for line in open(inputfilename)]

temperature = np.zeros(shape=(0, 0))
density = np.zeros(shape=(0, 0))

nt = 0  # get index of last unique variable before linebreak in input file

# get all unique temperatures
for i in range(1, len(aneosfile)):
    try:
        temperature = np.append(temperature, reformat(aneosfile[i][1]))  # indexing is row, column
    except IndexError:
        nt = i - 1
        break

# get all unique densities
for i in range(1, len(aneosfile), nt + 1):
    density = np.append(density, reformat(aneosfile[i][0]))  # indexing is row, column

# get number of unique densities
nr = len(density)  # density grid number

# create a matrix of #density X #variable
energy = np.zeros(shape=(nr, nt)) #J/kg
pressure = np.zeros(shape=(nr, nt)) #Pa
soundspeed = np.zeros(shape=(nr, nt)) #m/s
entropy = np.zeros(shape=(nr, nt)) #J/kg/K

"""
Recall that nr is the number of unique densities and nt is the number of unique samples of the variable in question.
"""
i = 1
for m in range(0, nr):
    for n in range(0, nt):
        try:
            # assign variable to row and column of corresponding matrix
            energy[m][n] = reformat(aneosfile[i][2])
            pressure[m][n] = reformat(aneosfile[i][3])
            soundspeed[m][n] = reformat(aneosfile[i][4])
            entropy[m][n] = reformat(aneosfile[i][5])

        except IndexError:  # skipping a line, I think for when you get to blank lines?
            i = i + 1
            energy[m][n] = reformat(aneosfile[i][2])
            pressure[m][n] = reformat(aneosfile[i][3])
            soundspeed[m][n] = reformat(aneosfile[i][4])
            entropy[m][n] = reformat(aneosfile[i][5])
        i = i + 1  # increase row number in input file


# Taking the min and max internal energy from the original ANEOS data
umin = np.min(energy)
umax = np.max(energy)

delta = (umax / umin)**(1.0 / (nu - 1))  # create exponential change for new grid

new_energy = np.zeros(shape=(0, 0))  # create a new energy array

"""
Create a new grid of energies that change exponentially away from the origin.
"""
for m in range(0, nu):
    new_energy = np.append(new_energy,umin*delta**m)  # exponential grid

# plt.rcParams["figure.figsize"] = [16, 9]
# plt.rcParams.update({'font.size': 16})
# fig_energy = plt.figure()
# ax_energy = fig_energy.add_subplot(111)
# ax_energy.plot(list(range(0, nu)), new_energy * 1e-6, linewidth=2.0, color='black')
# ax_energy.set_title('Modified Internal Energy vs. Grid Spacing')
# ax_energy.set_xlabel("Grid Number")
# ax_energy.set_ylabel("Modified Internal Energy (MJ/kg)")
# ax_energy.grid()

new_temperature = np.zeros(shape=(nr, nu))
new_pressure = np.zeros(shape=(nr, nu))
new_soundspeed = np.zeros(shape=(nr, nu))
new_entropy = np.zeros(shape=(nr, nu))


# 1D interpolation & extrapolation (linear)
for m in range(0, nu):
    # approximate temperature given internal energy
    f_temperature = interpolate.interp1d(energy[m, :], temperature, kind='linear', fill_value='extrapolate')
    new_temperature[m] = f_temperature(new_energy)

    # approximate pressure given temperature
    f_pressure = interpolate.interp1d(temperature, pressure[m, :], kind='linear', fill_value='extrapolate')
    new_pressure[m] = f_pressure(new_temperature[m][:])

    # approximate sound speed given temperature
    f_soundspeed = interpolate.interp1d(temperature, soundspeed[m, :], kind='linear', fill_value='extrapolate')
    new_soundspeed[m] = f_soundspeed(new_temperature[m][:])

    # approximate entropy given temperature
    f_entropy = interpolate.interp1d(temperature, entropy[m, :], kind='linear', fill_value='extrapolate')
    new_entropy[m] = f_entropy(new_temperature[m][:])



plt.rcParams["figure.figsize"] = [16, 9]
plt.rcParams.update({'font.size': 16})

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
for m in range(0, nr, int(nr/6)):
    d = density[m]
    ax1.semilogy(temperature, energy[m], label="{} kg/m3".format(float(round(d, 2))))
ax1.set_xlabel("Temperature (K)")
ax1.set_ylabel("Energy (J/Kg)")
ax1.set_title("Temperature vs. Energy For Selected Variable Densities in Database")
ax1.legend(loc='lower right')
ax1.grid()

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.semilogy([i for i in density], new_energy, linewidth=2, color='black')
ax2.set_xlabel("Density (kg/m3)")
ax2.set_ylabel("Modified Internal Energy (J/kg)")
ax2.set_title("Modified Internal Energy vs Density Grid")
ax2.grid()

plt.show()