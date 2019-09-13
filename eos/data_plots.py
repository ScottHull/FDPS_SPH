# This is a python script that converts u(rho, T), P(rho, T), Cs(rho,T), S(rho, T)
# to T(rho, u), P(rho, u), Cs(rho, u), S(rho, u), which is more useful for SPH calculaitons
#

import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.interpolate import interp1d
from scipy import interpolate

#----- A user has to change these three parameters  ----------------

inputfilename = "granite.table.txt"  # input ANEOS file. This follows the format from iSALE
outputfilename = "granite.rho_u.txt"  # output ANEOS file
nu = 120  # number of the grid for the internal energy (exponential)

#-------------------------------------------------------------------

# This function is to correct the original ANEOS format that does not include "E"
# This seems to  occur when the exponent reaches -101
def reformat(number):
    if number.find('E') == -1:
        exponent = "-101"
        mantissa  = number.split(exponent)
        return float(mantissa[0])*10**float(exponent)
    else:
        mantissa, exponent= number.split('E')

    return float(mantissa)*10**float(exponent)


aneosfile= [line.split() for line in open(inputfilename)]

temperature=np.zeros(shape=(0,0))
density=np.zeros(shape=(0,0))

for i in range(1, len(aneosfile)):
    try:
        temperature=np.append(temperature,reformat(aneosfile[i][1]))
    except IndexError:
        nt=i-1
        break

for i in range(1,len(aneosfile), nt+1):
    density=np.append(density,reformat(aneosfile[i][0]))


nr=len(density) #density grid number

energy=np.zeros(shape=(nr,nt)) #J/kg
pressure=np.zeros(shape=(nr,nt)) #Pa
soundspeed=np.zeros(shape=(nr,nt)) #m/s
entropy=np.zeros(shape=(nr,nt)) #J/kg/K

i = 1
for m in range(0,nr):
    for n in range(0,nt):
        try:
            energy[m][n] = reformat(aneosfile[i][2])
            pressure[m][n] = reformat(aneosfile[i][3])
            soundspeed[m][n] = reformat(aneosfile[i][4])
            entropy[m][n] = reformat(aneosfile[i][5])

        except IndexError: #skipping a line
            i = i+1
            energy[m][n] = reformat(aneosfile[i][2])
            pressure[m][n] = reformat(aneosfile[i][3])
            soundspeed[m][n] = reformat(aneosfile[i][4])
            entropy[m][n] = reformat(aneosfile[i][5])
        i = i+1


# Taking the min and max internal energy from the original ANEOS data
umin = np.min(energy)
umax = np.max(energy)

delta = (umax / umin)**(1.0 / (nu - 1))

new_energy = np.zeros(shape=(0, 0))
for m in range(0, nu):
    new_energy = np.append(new_energy, umin*delta**m) #exponential grid

new_temperature = np.zeros(shape=(nr, nu))
new_pressure = np.zeros(shape=(nr, nu))
new_soundspeed = np.zeros(shape=(nr, nu))
new_entropy = np.zeros(shape=(nr, nu))


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


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
for m in range(0, int(nu/6)):
    d = density[m]
    ax1.plot(temperature, energy[m] * 1e-6, label="{} kg/m3".format(d))
ax1.set_xlabel("Temperature (K)")
ax1.set_ylabel("Energy (MJ/kg)")
ax1.grid()

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
for m in range(0, int(nu/6)):
    d = density[m]
    ax2.plot3d(new_energy, [d for i in new_energy])

plt.show()
