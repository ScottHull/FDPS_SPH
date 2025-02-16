# This is a python script that converts u(rho, T), P(rho, T), Cs(rho,T), S(rho, T)
# to T(rho, u), P(rho, u), Cs(rho, u), S(rho, u), which is more useful for SPH calculaitons
#

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import sys
from scipy.interpolate import interp1d
from scipy import interpolate

# ----- A user has to change these three parameters  ----------------

inputfilename = "granite.table.txt"  # input ANEOS file. This follows the format from iSALE
outputfilename = "granite.rho_u.txt"  # output ANEOS file
nu = 120  # number of the grid for the internal energy (exponential)

#-------------------------------------------------------------------

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

plt.rcParams["figure.figsize"] = [16, 9]
plt.rcParams.update({'font.size': 16})
fig_energy = plt.figure()
ax_energy = fig_energy.add_subplot(111)
# ax_energy.plot(list(range(0, nu)), new_energy * 1e-6, linewidth=2.0, color='black')
ax_energy.plot(density, new_energy * 1e-6, linewidth=2.0, color="black")
ax_energy.set_title('Modified Internal Energy vs. Density')
ax_energy.set_xlabel("Density (kg/m3)")
ax_energy.set_ylabel("Modified Internal Energy (MJ/kg)")
ax_energy.grid()

new_temperature = np.zeros(shape=(nr, nu))
new_pressure = np.zeros(shape=(nr, nu))
new_soundspeed = np.zeros(shape=(nr, nu))
new_entropy = np.zeros(shape=(nr, nu))



# 1D interpolation & extrapolation (linear)
for m in range(0, nu):
    # approximate temperature given internal energy
    f_temperature = interpolate.interp1d(energy[m], temperature, kind='linear', fill_value='extrapolate')
    new_temperature[m] = f_temperature(new_energy)

    # approximate pressure given temperature
    f_pressure = interpolate.interp1d(energy[m], pressure[m], kind='linear', fill_value='extrapolate')
    new_pressure[m] = f_pressure(new_energy)

    # approximate sound speed given temperature
    f_soundspeed = interpolate.interp1d(energy[m], soundspeed[m], kind='linear', fill_value='extrapolate')
    new_soundspeed[m] = f_soundspeed(new_energy)

    # approximate entropy given temperature
    f_entropy = interpolate.interp1d(energy[m], entropy[m], kind='linear', fill_value='extrapolate')
    new_entropy[m] = f_entropy(new_energy)


# producing a few output images to make sure that this fitting is doing an okay job
# plt.rcParams["figure.figsize"] = [16, 9]
# for m in range(0, nr, int(nr/6)):
#
#     ax = [0, 0, 0, 0, 0]
#
#     fig = plt.figure(figsize=(10, 6.128))
#
#     ax[0] = fig.add_subplot(221)
#     ax[1] = fig.add_subplot(222)
#     ax[2] = fig.add_subplot(223)
#     ax[3] = fig.add_subplot(224)
#
#     ax[0].semilogy(energy[m] * 1e-6, temperature * 1e-3, '--', label="original ANEOS")
#     ax[0].semilogy(new_energy * 1e-6, new_temperature[m] * 1e-3, '-.', label="modified")
#
#
#     ax[1].semilogy(energy[m] * 1e-6, soundspeed[m] * 1e-3,'--')
#     ax[1].semilogy(new_energy * 1e-6, new_soundspeed[m] * 1e-3, '-.')
#
#     ax[2].semilogy(energy[m] * 1e-6, entropy[m] * 1e-3,'--')
#     ax[2].semilogy(new_energy * 1e-6, new_entropy[m] * 1e-3,'--')
#
#     ax[3].semilogy(energy[m] * 1e-6, pressure[m] * 1e-3,'--')
#     ax[3].semilogy(new_energy * 1e-6, new_pressure[m] * 1e-3,'--')
#
#
#     ax[0].legend(frameon=False)
#
#     ax[0].set_xlabel('Energy (MJ/kg)', fontsize=10)
#     ax[1].set_xlabel('Energy (MJ/kg)', fontsize=10)
#     ax[2].set_xlabel('Energy (MJ/kg)', fontsize=10)
#     ax[3].set_xlabel('Energy (MJ/kg)', fontsize=10)
#
#     ax[0].set_ylabel('Temperature (K)', fontsize=10)
#     ax[1].set_ylabel('Sound Speed (m/s2)', fontsize=10)
#     ax[2].set_ylabel('Entropy (kJ/K/kg)', fontsize=10)
#     ax[3].set_ylabel('Pressure (MPa)', fontsize=10)
#
#     ax[0].grid()
#     ax[1].grid()
#     ax[2].grid()
#     ax[3].grid()
#
#     fig.suptitle("Density: %3.3f kg/m$^3$" %(density[m]))
#     # fig.savefig("Density" + str(m) + ".png")




# new_temperature_2d = np.zeros(shape=(nr, nu, nu))
# new_pressure_2d = np.zeros(shape=(nr, nu, nu))
# new_soundspeed_2d = np.zeros(shape=(nr, nu, nu))
# new_entropy_2d = np.zeros(shape=(nr, nu, nu))

# 1D interpolation & extrapolation (linear)
# for m in range(0, nu):
#
#     # approximate temperature given internal energy
#     f_temperature = interpolate.interp2d(new_energy, density, new_temperature, kind='linear', fill_value='extrapolate')
#     new_temperature_2d[m] = f_temperature(new_energy, density)
#
#     # approximate pressure given temperature
#     f_pressure = interpolate.interp2d(new_energy, density, new_pressure, kind='linear', fill_value='extrapolate')
#     new_pressure_2d[m] = f_pressure(new_energy, density)
#
#     # approximate sound speed given temperature
#     f_soundspeed = interpolate.interp2d(new_energy, density, new_soundspeed, kind='linear', fill_value='extrapolate')
#     new_soundspeed_2d[m] = f_soundspeed(new_energy, density)
#
#     # approximate entropy given temperature
#     f_entropy = interpolate.interp2d(new_energy, density, new_entropy, kind='linear', fill_value='extrapolate')
#     new_entropy_2d[m] = f_entropy(new_energy, density)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for m in range(0, nr, int(nr/1)):
    ax.plot_surface(new_energy, density, new_temperature, rstride=1, cstride=1,
                    cmap='viridis')
    ax.plot(new_energy, density, linestyle='--', color='green')
    # ax.plot3D(energy[m], [density[m] for i in energy[m]], temperature, linestyle='--', color='red')
    # ax.plot3D(new_energy, density, new_temperature[m], linestyle='--', color='purple')
# for m in range(0, nr):
#     ax.plot3D(new_energy, density, new_temperature[m], linestyle='--', color='purple')
# for m in range(0, nu):
#     ax.plot3D(energy[m], [density[m] for i in energy[m]], temperature)
ax.set_xlabel('Energy (J/kg)')
ax.set_ylabel('Density (kg)')
ax.set_zlabel('Temperature (K)')


plt.show()


























# # 1D interpolation & extrapolation (linear)
# for m in range(0, nu):
#     # approximate temperature given internal energy
#     f_temperature = interpolate.interp1d(energy[m, :], temperature, kind='linear', fill_value='extrapolate')
#     new_temperature[m] = f_temperature(new_energy)
#
#     # approximate pressure given temperature
#     f_pressure = interpolate.interp1d(temperature, pressure[m, :], kind='linear', fill_value='extrapolate')
#     new_pressure[m] = f_pressure(new_temperature[m][:])
#
#     # approximate sound speed given temperature
#     f_soundspeed = interpolate.interp1d(temperature, soundspeed[m, :], kind='linear', fill_value='extrapolate')
#     new_soundspeed[m] = f_soundspeed(new_temperature[m][:])
#
#     # approximate entropy given temperature
#     f_entropy = interpolate.interp1d(temperature, entropy[m, :], kind='linear', fill_value='extrapolate')
#     new_entropy[m] = f_entropy(new_temperature[m][:])
#
#
# # producing a few output images to make sure that this fitting is doing an okay job
# plt.rcParams["figure.figsize"] = [16, 9]
# for m in range(0, nr, int(nr/6)):
#
#     ax = [0, 0, 0, 0]
#
#     fig = plt.figure(figsize=(10, 6.128))
#
#     ax[0] = fig.add_subplot(221)
#     ax[1] = fig.add_subplot(222)
#     ax[2] = fig.add_subplot(223)
#     ax[3] = fig.add_subplot(224)
#
#     ax[0].semilogy(temperature * 1e-3, energy[m] * 1e-6, '--', label="original ANEOS")
#     ax[0].semilogy(new_temperature[m] * 1e-3, new_energy * 1e-6, '-.', label="modified")
#     ax[1].semilogy(temperature * 1e-3, pressure[m] * 1e-6,'--', new_temperature[m] * 1e-3, new_pressure[m] * 1e-6, '-.')
#     ax[2].semilogy(temperature * 1e-3, soundspeed[m] * 1e-3,'--', new_temperature[m] * 1e-3, new_soundspeed[m] * 1e-3, '-.')
#     ax[3].semilogy(temperature * 1e-3, entropy[m] * 1e-3,'--', new_temperature[m] * 1e-3, new_entropy[m] * 1e-3, '-.')
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
#     ax[0].grid()
#     ax[1].grid()
#     ax[2].grid()
#     ax[3].grid()
#
#     fig.suptitle("Density: %3.3f kg/m$^3$" %(density[m]))
#     # fig.savefig("Density" + str(m) + ".png")
#
# plt.show()

# for m in range(0, nr, int(nr/6)):
#
#     ax = [0, 0, 0, 0]
#
#     fig = plt.figure(figsize=(10, 6.128))
#
#     ax[0] = fig.add_subplot(221)
#     ax[1] = fig.add_subplot(222)
#     ax[2] = fig.add_subplot(223)
#     ax[3] = fig.add_subplot(224)
#
#     ax[0].semilogy(energy[m] * 1e-6, temperature * 1e-3, '--', label="original ANEOS")
#     ax[0].semilogy(new_energy * 1e-6, new_temperature[m] * 1e-3, '-.', label="modified")
#     ax[1].semilogy(energy[m] * 1e-6, pressure[m] * 1e-6,'--', new_energy * 1e-6, new_pressure[m] * 1e-6, '-.')
#     ax[2].plot(energy[m] * 1e-6, soundspeed[m] * 1e-3,'--', new_energy * 1e-6, new_soundspeed[m] * 1e-3, '-.')
#     ax[3].plot(energy[m] * 1e-6, entropy[m] * 1e-3,'--', new_energy * 1e-6, new_entropy[m] * 1e-3, '-.')
#
#     ax[0].legend(frameon=False)
#
#     ax[0].set_xlabel('Energy (MJ/kg)', fontsize=10)
#     ax[1].set_ylabel('Pressure (MPa)', fontsize=10)
#     ax[2].set_ylabel('Sound Speed (km/s)', fontsize=10)
#     ax[3].set_ylabel('Entropy (kJ/K/kg)', fontsize=10)
#     ax[2].set_xlabel('Energy (MJ/kg)', fontsize=10)
#     ax[3].set_xlabel('Energy (MJ/kg)', fontsize=10)
#
#     ax[0].grid()
#     ax[1].grid()
#     ax[2].grid()
#     ax[3].grid()
#
#     fig.suptitle("Density: %3.3f kg/m$^3$" %(density[m]))
#     # fig.savefig("Density" + str(m) + ".png")
#
# for m in range(0, nr, int(nr/6)):
#
#     fig2 = plt.figure()
#     ax2 = fig2.add_subplot(111)
#     ax2.plot(energy * 1e-6, new_energ * 1e-6)
#     ax2.set_title("Internal Energy vs. Modified Internal Energy")
#     ax2.set_xlabel("Internal Energy (mJ/kg)")
#     ax2.set_ylabel("Modified Internal Energy (mJ/kg)")
#     ax2.grid()

# h = open(outputfilename,'w')
# h.write("%i %i %s \n" % (nr, nu , ':Grid numbers for density and internal energy'))
# h.write('Density (km/m3), Internal energy (kJ/kg), Temperature (K), Sound speed (m/s), Entropy (J/K/kg) \n')
#
#
# for m in range(0,nr):
#     for n in range(0,nu):
#         h.write("%15.8E %15.8E %15.8E %15.8E %15.8E \n" % (density[m], new_energy[n], new_temperature[m][n], new_soundspeed[m][n], new_entropy[m][n]))
