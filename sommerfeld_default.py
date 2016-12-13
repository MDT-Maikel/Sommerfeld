#! /usr/bin/env python

import sys
import argparse
import mpmath


##################
# xsec functions #
##################

def get_xsec(process, rep, l, sommerfeld, m, v, alphas, alphasommerfeld):
	if process == 'sstoqq':
		if rep == 3:
			return xsec_sstoqq_3(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 6:
			return xsec_sstoqq_6(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 8:
			return xsec_sstoqq_8(l, sommerfeld, m, v, alphas, alphasommerfeld)
		return 0.0

	if process == 'sstogg':
		if rep == 3:
			return xsec_sstogg_3(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 6:
			return xsec_sstogg_6(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 8:
			return xsec_sstogg_8(l, sommerfeld, m, v, alphas, alphasommerfeld)
		return 0.0;

	if process == 'fftoqq':
		if rep == 3:
			return xsec_fftoqq_3(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 6:
			return xsec_fftoqq_6(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 8:
			return xsec_fftoqq_8(l, sommerfeld, m, v, alphas, alphasommerfeld)
		return 0.0

	if process == 'fftogg':
		if rep == 3:
			return xsec_fftogg_3(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 6:
			return xsec_fftogg_6(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 8:
			return xsec_fftogg_8(l, sommerfeld, m, v, alphas, alphasommerfeld)
		return 0.0

	if process == 'vvtoqq':
		if rep == 3:
			return xsec_vvtoqq_3(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 6:
			return xsec_vvtoqq_6(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 8:
			return xsec_vvtoqq_8(l, sommerfeld, m, v, alphas, alphasommerfeld)
		return 0.0

	if process == 'vvtogg':
		if rep == 3:
			return xsec_vvtogg_3(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 6:
			return xsec_vvtogg_6(l, sommerfeld, m, v, alphas, alphasommerfeld)
		if rep == 8:
			return xsec_vvtogg_8(l, sommerfeld, m, v, alphas, alphasommerfeld)
		return 0.0

	return 0.0
	

def xsec_sstoqq_3(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_S3S3_QQ_NS_L0, ID_S3S3_QQ_NS_L1, ID_S3S3_QQ_NS_L2, ID_S3S3_QQ_NS_L3, ID_S3S3_QQ_NS_L4]
	else:
		wave_list = [ID_S3S3_QQ_SO_L0, ID_S3S3_QQ_SO_L1, ID_S3S3_QQ_SO_L2, ID_S3S3_QQ_SO_L3, ID_S3S3_QQ_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_sstoqq_6(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_S6S6_QQ_NS_L0, ID_S6S6_QQ_NS_L1, ID_S6S6_QQ_NS_L2, ID_S6S6_QQ_NS_L3, ID_S6S6_QQ_NS_L4]
	else:
		wave_list = [ID_S6S6_QQ_SO_L0, ID_S6S6_QQ_SO_L1, ID_S6S6_QQ_SO_L2, ID_S6S6_QQ_SO_L3, ID_S6S6_QQ_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_sstoqq_8(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_S8S8_QQ_NS_L0, ID_S8S8_QQ_NS_L1, ID_S8S8_QQ_NS_L2, ID_S8S8_QQ_NS_L3, ID_S8S8_QQ_NS_L4]
	else:
		wave_list = [ID_S8S8_QQ_SO_L0, ID_S8S8_QQ_SO_L1, ID_S8S8_QQ_SO_L2, ID_S8S8_QQ_SO_L3, ID_S8S8_QQ_SO_L4]
	return sum(wave_list[:l + 1])
	

def xsec_sstogg_3(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_S3S3_GG_NS_L0, ID_S3S3_GG_NS_L1, ID_S3S3_GG_NS_L2, ID_S3S3_GG_NS_L3, ID_S3S3_GG_NS_L4]
	else:
		wave_list = [ID_S3S3_GG_SO_L0, ID_S3S3_GG_SO_L1, ID_S3S3_GG_SO_L2, ID_S3S3_GG_SO_L3, ID_S3S3_GG_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_sstogg_6(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_S6S6_GG_NS_L0, ID_S6S6_GG_NS_L1, ID_S6S6_GG_NS_L2, ID_S6S6_GG_NS_L3, ID_S6S6_GG_NS_L4]
	else:
		wave_list = [ID_S6S6_GG_SO_L0, ID_S6S6_GG_SO_L1, ID_S6S6_GG_SO_L2, ID_S6S6_GG_SO_L3, ID_S6S6_GG_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_sstogg_8(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_S8S8_GG_NS_L0, ID_S8S8_GG_NS_L1, ID_S8S8_GG_NS_L2, ID_S8S8_GG_NS_L3, ID_S8S8_GG_NS_L4]
	else:
		wave_list = [ID_S8S8_GG_SO_L0, ID_S8S8_GG_SO_L1, ID_S8S8_GG_SO_L2, ID_S8S8_GG_SO_L3, ID_S8S8_GG_SO_L4]
	return sum(wave_list[:l + 1])


def xsec_fftoqq_3(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_F3F3_QQ_NS_L0, ID_F3F3_QQ_NS_L1, ID_F3F3_QQ_NS_L2, ID_F3F3_QQ_NS_L3, ID_F3F3_QQ_NS_L4]
	else:
		wave_list = [ID_F3F3_QQ_SO_L0, ID_F3F3_QQ_SO_L1, ID_F3F3_QQ_SO_L2, ID_F3F3_QQ_SO_L3, ID_F3F3_QQ_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_fftoqq_6(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_F6F6_QQ_NS_L0, ID_F6F6_QQ_NS_L1, ID_F6F6_QQ_NS_L2, ID_F6F6_QQ_NS_L3, ID_F6F6_QQ_NS_L4]
	else:
		wave_list = [ID_F6F6_QQ_SO_L0, ID_F6F6_QQ_SO_L1, ID_F6F6_QQ_SO_L2, ID_F6F6_QQ_SO_L3, ID_F6F6_QQ_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_fftoqq_8(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_F8F8_QQ_NS_L0, ID_F8F8_QQ_NS_L1, ID_F8F8_QQ_NS_L2, ID_F8F8_QQ_NS_L3, ID_F8F8_QQ_NS_L4]
	else:
		wave_list = [ID_F8F8_QQ_SO_L0, ID_F8F8_QQ_SO_L1, ID_F8F8_QQ_SO_L2, ID_F8F8_QQ_SO_L3, ID_F8F8_QQ_SO_L4]
	return sum(wave_list[:l + 1])


def xsec_fftogg_3(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_F3F3_GG_NS_L0, ID_F3F3_GG_NS_L1, ID_F3F3_GG_NS_L2, ID_F3F3_GG_NS_L3, ID_F3F3_GG_NS_L4]
	else:
		wave_list = [ID_F3F3_GG_SO_L0, ID_F3F3_GG_SO_L1, ID_F3F3_GG_SO_L2, ID_F3F3_GG_SO_L3, ID_F3F3_GG_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_fftogg_6(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_F6F6_GG_NS_L0, ID_F6F6_GG_NS_L1, ID_F6F6_GG_NS_L2, ID_F6F6_GG_NS_L3, ID_F6F6_GG_NS_L4]
	else:
		wave_list = [ID_F6F6_GG_SO_L0, ID_F6F6_GG_SO_L1, ID_F6F6_GG_SO_L2, ID_F6F6_GG_SO_L3, ID_F6F6_GG_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_fftogg_8(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_F8F8_GG_NS_L0, ID_F8F8_GG_NS_L1, ID_F8F8_GG_NS_L2, ID_F8F8_GG_NS_L3, ID_F8F8_GG_NS_L4]
	else:
		wave_list = [ID_F8F8_GG_SO_L0, ID_F8F8_GG_SO_L1, ID_F8F8_GG_SO_L2, ID_F8F8_GG_SO_L3, ID_F8F8_GG_SO_L4]
	return sum(wave_list[:l + 1])


def xsec_vvtoqq_3(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_V3V3_QQ_NS_L0, ID_V3V3_QQ_NS_L1, ID_V3V3_QQ_NS_L2, ID_V3V3_QQ_NS_L3, ID_V3V3_QQ_NS_L4]
	else:
		wave_list = [ID_V3V3_QQ_SO_L0, ID_V3V3_QQ_SO_L1, ID_V3V3_QQ_SO_L2, ID_V3V3_QQ_SO_L3, ID_V3V3_QQ_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_vvtoqq_6(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_V6V6_QQ_NS_L0, ID_V6V6_QQ_NS_L1, ID_V6V6_QQ_NS_L2, ID_V6V6_QQ_NS_L3, ID_V6V6_QQ_NS_L4]
	else:
		wave_list = [ID_V6V6_QQ_SO_L0, ID_V6V6_QQ_SO_L1, ID_V6V6_QQ_SO_L2, ID_V6V6_QQ_SO_L3, ID_V6V6_QQ_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_vvtoqq_8(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_V8V8_QQ_NS_L0, ID_V8V8_QQ_NS_L1, ID_V8V8_QQ_NS_L2, ID_V8V8_QQ_NS_L3, ID_V8V8_QQ_NS_L4]
	else:
		wave_list = [ID_V8V8_QQ_SO_L0, ID_V8V8_QQ_SO_L1, ID_V8V8_QQ_SO_L2, ID_V8V8_QQ_SO_L3, ID_V8V8_QQ_SO_L4]
	return sum(wave_list[:l + 1])


def xsec_vvtogg_3(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_V3V3_GG_NS_L0, ID_V3V3_GG_NS_L1, ID_V3V3_GG_NS_L2, ID_V3V3_GG_NS_L3, ID_V3V3_GG_NS_L4]
	else:
		wave_list = [ID_V3V3_GG_SO_L0, ID_V3V3_GG_SO_L1, ID_V3V3_GG_SO_L2, ID_V3V3_GG_SO_L3, ID_V3V3_GG_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_vvtogg_6(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_V6V6_GG_NS_L0, ID_V6V6_GG_NS_L1, ID_V6V6_GG_NS_L2, ID_V6V6_GG_NS_L3, ID_V6V6_GG_NS_L4]
	else:
		wave_list = [ID_V6V6_GG_SO_L0, ID_V6V6_GG_SO_L1, ID_V6V6_GG_SO_L2, ID_V6V6_GG_SO_L3, ID_V6V6_GG_SO_L4]
	return sum(wave_list[:l + 1])

def xsec_vvtogg_8(l, sommerfeld, m, v, alphas, alpha_sommerfeld):
	wave_list = []
	if not sommerfeld:
		wave_list = [ID_V8V8_GG_NS_L0, ID_V8V8_GG_NS_L1, ID_V8V8_GG_NS_L2, ID_V8V8_GG_NS_L3, ID_V8V8_GG_NS_L4]
	else:
		wave_list = [ID_V8V8_GG_SO_L0, ID_V8V8_GG_SO_L1, ID_V8V8_GG_SO_L2, ID_V8V8_GG_SO_L3, ID_V8V8_GG_SO_L4]
	return sum(wave_list[:l + 1])


###############
# main script #
###############

# argument parser
parser = argparse.ArgumentParser(description='Returns the cross section (in 1/GeV^2) of the annihilation process for a given partial wave.')
parser.add_argument('-p', '--process', action='store', required=True, help='the annihilation process (sstoqq, sstogg, fftoqq, fftogg, vvtoqq or vvtogg)')
parser.add_argument('-r', '--rep', action='store', required=True, help='the color representation of the annihilating particle (3, 6 or 8)')
parser.add_argument('-v', '--vars', nargs=4, required=True, help='the variables the cross section depends on: m, v, alpha_s, alpha_sommerfeld')
parser.add_argument('-s', '--sommerfeld', action='store_true', help='add sommerfeld corrections (default: off)')
parser.add_argument('-l', '--lwave', action='store', default=2, help='corrections up to the lth partial wave (default l = 2)')
parser.add_argument('-e', '--extended', action='store_true', help='Extended logging of the results')
args = parser.parse_args()

# determine the process and representation
process = args.process
if not process in ['sstoqq', 'sstogg', 'fftoqq', 'fftogg', 'vvtoqq', 'vvtogg']:
	print "Process " + process + " is not known, must be sstoqq, sstogg, fftoqq, fftogg, vvtoqq or vvtogg."
	sys.exit(2)
rep = int(args.rep)
if not rep in [3, 6, 8]:
	print "Color representation " + str(rep) + " is not valid, must be 3, 6, or 8."
	sys.exit(2)

# determine variables
m = mpmath.mpf(args.vars[0])
v = mpmath.mpf(args.vars[1])
alphas = mpmath.mpf(args.vars[2])
alphasommerfeld = mpmath.mpf(args.vars[3])

# determine partial wave and sommerfeld
l = int(args.lwave);
sommerfeld = args.sommerfeld

# determine the cross section
xsec = get_xsec(process, rep, l, sommerfeld, m, v, alphas, alphasommerfeld);

# print the result
if args.extended:
	print "Annihilation cross section for " + process + " with color representation " + str(rep)
	print "and m = " + str(m) + ", v = " + str(v) + ", alpha_s = " + str(alphas) + ", alpha_sommerfeld = " + str(alphasommerfeld)
	print "and for l = " + str(l) + " and sommerfeld " + str(sommerfeld) + " equals:"
	print "\t" + mpmath.nstr(xsec, 15)
else:
	print mpmath.nstr(xsec, 15)

