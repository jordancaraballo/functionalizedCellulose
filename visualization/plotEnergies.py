#!/bin/python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------------------------------
# Script designed to plot energies output (specially potential and cinetic) from NAMD and LAMMPS simulations.
# The script is designed to create a simple plot. There is nothing complicated here. It works to facilitate and
# speed up other processes.
# Author: Jordan A Caraballo-Vega

import matplotlib.pyplot as plt
import numpy as np
import sys

# Legend of namd2 output file
# ETITLE: TS, BOND, ANGLE, DIHED, IMPRP, ELECT, VDW, BOUNDARY, MISC, KINETIC, TOTAL,
# TEMP, POTENTIAL, TOTAL3, TEMPAVG, PRESSURE, GPRESSURE, VOLUME, PRESSAVG, GPRESSAVG

file_name  = open(sys.argv[1], 'r')
title_name = open(sys.argv[2], 'r')
energies   = [line.strip().split() for line in open(file_name,'r').readlines() ][9:]

if not title_name:
    title_name = "Add Simulation Description Here"

t = list()
s = list()
for i in energies:
    t.append(i[1])
    s.append(i[13])

plt.plot(t, s)
plt.xlabel('Time (ts)')
plt.ylabel('Potential (J)')
plt.suptitle('Potential as Function of Time', fontsize=16)
plt.title(title_name)
plt.grid(True)
plt.savefig(title_name+".png")
plt.show()
