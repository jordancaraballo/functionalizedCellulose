# -*- coding: utf-8 -*-
# pruebaDCD.py
# tests DCDReader from Wolffia's library
# git pull wolffia first
# use: python pruebaDCD.py <path to dcd file>
#
# José O. Sotero Esteva
# UPRH Humacao, February 6, 2018
#


if __name__ == '__main__':
	import sys, os
	sys.path.append('/home/jse/inv/Cuchifritos/git/wolffia/')

	from lib.io.dcd import DCDReader

	print sys.argv[1]
	dcd = DCDReader(sys.argv[1])

	positions = [12,20,30]
	print "Coordinates in positions ", positions

	for frame in dcd:
		print [(frame[0][p], frame[1][p], frame[2][p]) for p in positions]
