#! /usr/bin/env python

# python modules
import sys
import numpy
import os
import subprocess
import math

# create output file
outputfile = open("rd_mass.txt", "w")
outputfile.write("mass_dm mass_x delta omega omega(sommerfeld) omega(bsf) omega(sommerfeld + bsf)\n")

# initialize value for delta
delta = 0.0

# loop over dark matter mass
for mdm in numpy.arange(1.0, 6001.0 + 0.1, 25.0):
	print "mdm = ", mdm

	# update parameters
	mx = mdm * (1.0 + delta)

	# create micromegas param file
	params = "MDM "+ str(mdm) + "\nMX " + str(mx)
	with open("input_micromegas.par", 'w') as param_file:
		param_file.write(params)

	# run main micromegas
	output = subprocess.check_output("./main input_micromegas.par", shell=True)
	omega = float([line for line in output.split("\n") if "omega_h^2 = " in line].pop().split().pop())

	# run main micromegas with sommerfeld corrections
	output_sommerfeld = subprocess.check_output("./main input_micromegas.par on", shell=True)
	omega_sommerfeld = float([line for line in output_sommerfeld.split("\n") if "omega_h^2 = " in line].pop().split().pop())


	# run main micromegas with bound state corrections
	output_bsf = subprocess.check_output("./main input_micromegas.par off on", shell=True)
	omega_bsf = float([line for line in output_bsf.split("\n") if "omega_h^2 = " in line].pop().split().pop())

	# run main micromegas with sommerfeld and bound state corrections
	output_sommerfeld_bsf = subprocess.check_output("./main input_micromegas.par on on", shell=True)
	omega_sommerfeld_bsf = float([line for line in output_sommerfeld_bsf.split("\n") if "omega_h^2 = " in line].pop().split().pop())

	# write to output file
	outputfile.write("%.4f %.4f %.4f %.4e %.4e %.4e %.4e\n" % (mdm, mx, delta, omega, omega_sommerfeld, omega_bsf, omega_sommerfeld_bsf))
