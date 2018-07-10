#!/usr/bin/env python3

# Data structure: 
# 	distn
# 		mean
# 		sd
# 		ptile
# 	param: value pairs
# 	fitness: distn
# 	performance: distn
# 	genotype: distn list of loci + corr
# 	stoch: distn list of loci + corr + xcorr
	 
import os
import sys
import re
from inspect import currentframe, getframeinfo

def main():
	[infile, outfile] = read_args_open_files()
	read_runs(infile, outfile)
	close_files(infile, outfile)

def read_args_open_files():
	if len(sys.argv) != 2:
		print("\n\tUsage: extractData.py dataFile\n")
		exit()
	if not os.path.isfile(sys.argv[1]):
		print("\n\tFile {} not found".format(sys.argv[1]))
		exit()
	infile = open(sys.argv[1], 'r')
	outfile = open(sys.argv[1] + ".mma", 'w')
	outfile.write("Dataset[{\n")
	return [infile, outfile]
	
def close_files(infile, outfile):
	outfile.write("\n}]\n")
	infile.close()
	outfile.close()
	
def read_runs(infile, outfile):
	firstrun = True
	param = {}
	fdistn = {}
	pdistn = {}
	for line in infile:
		fields = line.split()
		# param group
		if not line.isspace() and fields[0] == "Run":
			print("{} {} {}".format(fields[0], fields[1], fields[2]))
			group = 1
			if not firstrun:
				print_data(param, fdistn, pdistn, outfile, False)
				param = {}
				# reset all dicts
			firstrun = False
		if group == 1:
			if not line.isspace():
				param[fields[0]] = fields[2]
			else:
				group = 2
				start = False
				ptile = []
				gdistn = []
				gtile = []
				for i in range(int(param["loci"])):
					gdistn.append({})
					gtile.append([1])
			continue
		# fitness distn group
		if group == 2:
			[group,start,ptile] = set_distn(fdistn,ptile,line,start,group,1)
			continue
		# performance distn group
		if group == 3:
			[group,start,ptile] = set_distn(pdistn,ptile,line,start,group,1)
			continue
		# genotype distn group
		if group == 4:
			for i in range(int(param["loci"])):
				[group,start,gtile[i]] = set_distn(gdistn[i],gtile[i],line,start,group,i)
			if group > 4: group = 5
			continue
			
	print_data(param, fdistn, pdistn, outfile, True)

def set_distn(distn, ptile, line, start, group, field_num):
	fields = line.split()
	if not line.isspace() and (fields[0] == "Mean" or fields[0] == "SD"):
		distn[fields[0]] = fields[field_num]
		start = True
	elif not line.isspace() and start:
		ptile.append([float(fields[0]),
			str(float(fields[1])).replace("e", "*10^")])
	elif line.isspace and start:
		distn["ptile"] = ptile
		group += 1
		start = False
		ptile = []
	return [group, start, ptile]

def print_param(param, outfile):
	outfile.write("<|\"param\" -> <|")
	first = True
	for k in param:
		if not first:
			outfile.write(",")
		first = False
		# change e -> E because MMA only recognizes E for sci notation
		outfile.write("\"{}\" -> {}".format(k, param[k].replace("e", "*10^")))
	outfile.write("|>|>,")

def print_distn(distn, name, outfile):
	outfile.write(("<|\"{}\" -> <| \"mean\" -> {}, \"sd\" -> {}" +
			", \"ptile\" -> {}|>|>").format(name, (distn["Mean"]).replace("e", "*10^"), 
			(distn["SD"]).replace("e", "*10^"), 
			str(distn["ptile"]).replace("[", "{").replace("]", "}").replace("'","")))

def print_data(param, fdistn, pdistn, outfile, last):
	outfile.write("<|")
	print_param(param, outfile)
	print_distn(fdistn, "fdistn", outfile)
	outfile.write(",")
	print_distn(pdistn, "pdistn", outfile)
	outfile.write("|>")
	if not last:
		outfile.write(",")
	
main()


# Dataset[{
# 	<|"run" -> 1, "fitD" -> <|"mean" -> 2, "sd" -> 3, 
#      "ptile" -> {{0, 1}, {1, 2}}|>|>,
#      <|"run" -> 1, "fitD" -> <|"mean" -> 2, "sd" -> 3, 
#      "ptile" -> {{0, 1}, {1, 2}}|>|>
# }]
