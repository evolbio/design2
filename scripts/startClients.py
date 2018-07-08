#!/usr/bin/env python3

from subprocess import check_output
import os
import argparse

job = "/Users/steve/sim/plasticity/02.sensitivity/scripts/jobRepeat.bash"
path = '/Users/steve/sim/plasticity/02.sensitivity/'
command = "sensitivity_client"
server = "fisher:50051"
allhosts = ["fisher", "alice2", "rex", "ficus", "batty", "localhost"]

parser = argparse.ArgumentParser(description='Start sim clients on list of hosts.',
	usage = 'startClients.py [-h] [-a] [-n num] [-f file] [hosts ...]',
	epilog = 'startClients.py [-a] [-n num] [-f file] [hosts], where no args returns help; \
	-a adds default hosts to list, and hosts is a list of hosts; if hosts and no -a, \
	then use only hosts')
parser.add_argument("-a", "--all", action='store_true', help="all hosts")
parser.add_argument("-n", "--num", type=int, default=1, help="num proc per host")
parser.add_argument("-f", "--file", help="read hosts from file, overrides others")
parser.add_argument("hosts", nargs=argparse.REMAINDER, help='list of hosts')
args = parser.parse_args()

if not args.all and not args.file and args.hosts == []:
	parser.print_help()
	print()
	exit()
	
if not args.file:
	if not args.all:
		hosts = args.hosts
	else:
		hosts = allhosts + args.hosts
else:
	hosts = []
	chosts = []
	if not os.path.isfile(args.file):
		print("\n\tFile {} not found".format(args.file))
		exit()
	f = open(args.file, 'r')
	for line in f:
		if line[0] != '#' and not line.isspace():
			[host, num] = line.split()
			hosts.append([host,num])
			
if hosts == []:
	print("\n\tNo hosts in hostfile")

num = args.num
for h in hosts:
	if args.file:
		num = h[1]
		h = h[0]
		chosts.append(h)
	print("Starting {} jobs on {}".format(num, h))
	check_output(['ssh', h, job, str(num), path + command, server])

if args.file: hosts = chosts
print("\nChecking status ...\n")
hosts = list(set(hosts))	# removes duplicates from list
output = check_output([path + 'scripts/statusClients.py'] + hosts)
print(output.decode("utf-8"), end='')