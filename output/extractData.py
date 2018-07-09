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
	for line in infile:
		fields = line.split()
		if not line.isspace() and fields[0] == "Run":
			print("{} {} {}".format(fields[0], fields[1], fields[2]))
			group = 1
			if not firstrun:
				print_data(param, fdistn, outfile, False)
				param = {}
				# reset all dicts
			firstrun = False
		if group == 1:
			if not line.isspace():
				param[fields[0]] = fields[2]
			else:
				group = 2
				start = False
				fdistn = {}
				distn = []
			continue
		if group == 2:
			if not line.isspace() and (fields[0] == "Mean" or fields[0] == "SD"):
				fdistn[fields[0]] = fields[1]
				start = True
			elif not line.isspace() and start:
				distn.append([float(fields[0]),
					str(float(fields[1])).replace("e", "*10^")])
			elif line.isspace and start:
				fdistn["distn"] = distn
				group = 3
				start = False
				pdistn = {}
				distn = []
			
	print_data(param, fdistn, outfile, True)
	
#def print_distn(distn, name, outfile):
	
	
def print_data(param, fdistn, outfile, last):
	outfile.write("<|")
	outfile.write("<|\"param\" -> <|")
	first = True
	for k in param:
		if not first:
			outfile.write(",")
		first = False
		# change e -> E because MMA only recognizes E for sci notation
		outfile.write("\"{}\" -> {}".format(k, param[k].replace("e", "*10^")))
	outfile.write("|>|>,")
	outfile.write(("<|\"fdistn\" -> <| \"mean\" -> {}, \"sd\" -> {}" +
			", \"distn\" -> {}|>|>|>").format((fdistn["Mean"]).replace("e", "*10^"), 
			(fdistn["SD"]).replace("e", "*10^"), 
			str(fdistn["distn"]).replace("[", "{").replace("]", "}").replace("'","")))
	if not last:
		outfile.write(",")
	
main()


# newrun = False
# firstrun = True
# for line in infile:
# 	fields = line.split()
# 	if not line.isspace() and fields[0] is "Run":
# 		param = None
# 		fitness = None
# 		perform = None
# 		genotype = None
# 		stoch = None
# 		if not firstrun:
# 			outfile.write(",\n")
# 		else:
# 			firstrun = False
# 		setnum = 0
# 		run = fields[2]
# 	if startparam:
# 		if line.isspace(): 
# 			setnum = False
# 		elif fields[0] not "recT":
# 			param.append()


# Dataset[{
# 	<|"run" -> 1, "fitD" -> <|"mean" -> 2, "sd" -> 3, 
#      "ptile" -> {{0, 1}, {1, 2}}|>|>,
#      <|"run" -> 1, "fitD" -> <|"mean" -> 2, "sd" -> 3, 
#      "ptile" -> {{0, 1}, {1, 2}}|>|>
# }]
