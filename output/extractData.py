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
import bz2
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
	infile_name = sys.argv[1]
	if "bz2" in infile_name:
		infile = bz2.open(infile_name, 'rt')
	else:
		infile = open(infile_name, 'rt')
	infile_name = infile_name.replace('.bz2', '')
	outfile = bz2.open(infile_name + ".mma.bz2", 'wt')
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
				print_data(param,fdistn,pdistn,gdistn,gcorr,sdistn,\
					scorr,sgcorr,outfile,stoch,False)
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
				sdistn = []
				gtile = []
				gcorr = []
				scorr = []
				sgcorr = []
				loci = int(param["loci"])
				stoch = True if float(param["stochWt"]) > 1e-6 else False
				for i in range(loci):
					gdistn.append({})
					sdistn.append({})
					gtile.append([])
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
			for i in range(loci):
				[group,start,gtile[i]] = \
					set_distn(gdistn[i],gtile[i],line,start,group,i+1,loci)
			if group > 4: group = 5
			continue
		# genotype corr group
		if group == 5:
			[group,start] = set_corr(gcorr,line,start,group,loci)
			continue
		# stoch genotype distn group
		if group == 6:
			for i in range(loci):
				[group,start,gtile[i]] = \
					set_distn(sdistn[i],gtile[i],line,start,group,i+1,loci)
			if group > 6: group = 7
			continue
		# stochastic corr group
		if group == 7:
			[group,start] = set_corr(scorr,line,start,group,loci)
			continue
		# stochastic X genotype corr group
		if group == 8:
			[group,start] = set_corr(sgcorr,line,start,group,loci)
			continue
			
	print_data(param,fdistn,pdistn,gdistn,gcorr,sdistn,scorr,sgcorr,outfile,stoch,True)

def set_distn(distn, ptile, line, start, group, field_num, field_max=1):
	fields = line.split()
	if not line.isspace() and (fields[0] == "Mean" or fields[0] == "SD"):
		distn[fields[0]] = fields[field_num]
		if field_num == field_max: start = True
	elif not line.isspace() and start:
		ptile.append([float(fields[0]),
			str(float(fields[field_num])).replace("e", "*10^")])
	elif line.isspace and start:
		distn["ptile"] = ptile
		group += 1
		if field_num == field_max: start = False
		ptile = []
	return [group, start, ptile]
	
def set_corr(corr, line, start, group, loci):
	fields = line.split()
	if not line.isspace() and ("." in fields[1]):
		corr.append(fields[1:(loci+1)])
		return [group,True]
	elif line.isspace and start:
		return [group+1,False]
	return [group,start]

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
			
def print_gdistn(gdistn, name, outfile):
	outfile.write("<|\"{}\" -> ".format(name) + "<|")
	for i in range(len(gdistn)):
		print_distn(gdistn[i], "g" + str(i), outfile)
		if i != (len(gdistn) - 1): outfile.write(",")
	outfile.write("|>|>")

def print_corr(corr, name, outfile):
	outfile.write("<|\"{}\" -> ".format(name))
	outfile.write("{}".format(corr).replace("[", "{").replace("]", "}").replace("'",""))
	outfile.write("|>")

def print_data(param,fdistn,pdistn,gdistn,gcorr,sdistn,scorr,sgcorr,outfile,stoch,last):
	outfile.write("<|")
	print_param(param, outfile)
	print_distn(fdistn, "fdistn", outfile)
	outfile.write(",")
	print_distn(pdistn, "pdistn", outfile)
	outfile.write(",")
	print_gdistn(gdistn, "gdistn", outfile)
	outfile.write(",")
	print_corr(gcorr, "gcorr", outfile)
	if stoch:
		outfile.write(",")
		print_gdistn(sdistn, "sdistn", outfile)		
		outfile.write(",")
		print_corr(scorr, "scorr", outfile)		
		outfile.write(",")
		print_corr(sgcorr, "sgcorr", outfile)		
	outfile.write("|>")
	if not last:
		outfile.write(",")
	
main()





