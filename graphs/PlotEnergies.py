#!/bin/python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------------------------------
# Script designed to plot energies output (specially potential and cinetic) from NAMD and LAMMPS simulations.
# The script is designed to create a simple plot. There is nothing complicated here. It works to facilitate and
# speed up other processes.

import matplotlib.pyplot as plt
import numpy as np

# Legend of output file
#ETITLE: TS, BOND, ANGLE, DIHED, IMPRP, ELECT, VDW, BOUNDARY, MISC, KINETIC, TOTAL, TEMP, POTENTIAL, TOTAL3, TEMPAVG, PRESSURE, GPRESSURE, VOLUME, PRESSAVG, GPRESSAVG

file_name = "Pentane100_EnergyOutput.txt"
energies  = [line.strip().split() for line in open(file_name,'r').readlines() ][9:]       # list with psf file lines

t = list()
s = list()
for i in energies:
    t.append(i[1])
    s.append(i[13])

plt.plot(t, s)

plt.xlabel('Time (ts)')
plt.ylabel('Potential (J)')
plt.suptitle('Potential as Function of Time', fontsize=16)
plt.title('100% 5-atoms spacer arms')
plt.grid(True)
plt.savefig("potential.png")
plt.show()
