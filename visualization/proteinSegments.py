#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 08:48:30 2019

@author: jse

Script to fix protein segments on coordinate and pdb files. By fixing these segments
you will be able to visualize ribbons on vmd.
"""

import sys

aminos =  {'ALA': 11, 'ARG': 25, 'ASN': 15, 'ASP': 13, 'CYS': 11, 'GLN': 18, \
           'GLU': 16 , 'GLY': 15, 'HSD': 18, 'ILE': 20, 'LEU': 20, 'LYS': 25, 'OOR': 2, \
           'MET': 18, 'PHE': 21, 'PRO': 16, 'SER': 12, 'THR': 15, 'TRP': 25, 'TYR': 22, 'VAL': 17}

cont = 1
seg = "WAT"
i=0

f = open(sys.argv[1], 'r')
for linea in f:
    #print('"'+ linea[17:20]+'"')
    if linea[17:20] in aminos:
        if seg != linea[17:20]: 
            #print(seg + ': ' + str(i))
            cont += 1
            seg = linea[17:20]
            i = 1
        elif i == aminos[seg]:
            i = 1
            cont += 1
        print((linea[:23] + "{:3d}" + linea[26:]).format(cont), end='')
        i += 1
    else:
        print(linea, end='')
