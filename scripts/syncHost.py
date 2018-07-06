#!/usr/bin/env python3

from subprocess import check_output
import argparse

command = ['ps', 'x', '-o', 'pid,%cpu,%mem,command', '|', 'grep']
allhosts = ["fisher", "alice2", "rex", "ficus", "batty", "localhost"]

parser = argparse.ArgumentParser(description='Show running clients on hosts.',
	epilog = 'syncHost.py [-n] [-a] [hosts], -a adds default hosts to list, and \
	hosts is the list of hosts, if no -a  and hosts, then use only hosts')
parser.add_argument("-s", "--silent", action='store_true', help="no rsync output")
parser.add_argument("-a", "--all", action='store_true', help="all hosts")
parser.add_argument("hosts", nargs=argparse.REMAINDER, help='list of hosts')
args = parser.parse_args()

if not args.all and args.hosts == []:
	parser.print_help()
	print()
	exit()
if not args.all:
	hosts = args.hosts
else:
	hosts = allhosts + args.hosts

for h in hosts:
	command = ['rsync', '-aNHAXxvz', '--delete', '/Users/steve/sim/simlib',
				'/Users/steve/sim/grpcControl', h + ':sim/.']
	output = check_output(command)
	command = ['rsync', '-aNHAXxvz', '--delete',
				'/Users/steve/sim/plasticity/02.sensitivity',
				'/Users/steve/sim/grpcControl', h + ':sim/plasticity/.']
	output += check_output(command)
	print(h + ':')
	if not args.silent:
		for line in output.splitlines():
			print("  " + line.decode("utf-8"))
