#!/usr/bin/env python3

from subprocess import check_output
import argparse
import re

proc = 'sensitivity_client'
command = ['ps', 'x', '-o', 'pid,%cpu,%mem,command', '|', 'grep']
allhosts = ["fisher", "alice2", "rex", "ficus", "batty", "localhost"]

parser = argparse.ArgumentParser(description='Show running clients on hosts.',
	epilog = 'statusClients.py [-p PROC] [-a] [hosts], where no args is all default\
	 hosts, -a adds default hosts to list, and hosts is the list of hosts,\
	  if no -a  and hosts, then use only hosts')
parser.add_argument("-p", "--proc", help="process name to match")
parser.add_argument("-a", "--all", action='store_true', help="all hosts")
parser.add_argument("hosts", nargs=argparse.REMAINDER, help='list of hosts')
args = parser.parse_args()

if (args.all or args.hosts == []):
	hosts = allhosts + args.hosts
else:
	hosts = args.hosts
if (args.proc != None):
	proc = args.proc

#print ("Running:", ' '.join(['ssh', 'host'] + command + [proc]))

for h in hosts:
	output = check_output(['ssh', h] + command + [proc])
	print(h + ':')
	for line in output.splitlines():
		line = line.decode("utf-8")
		if not 'statusClient' in line and not 'grep' in line:
			line = re.sub('/.*/', '', line)
			print("  " + line)
print()