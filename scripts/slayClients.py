#!/usr/bin/env python3

from subprocess import check_output
import argparse

command = "sensitivity_client"
allhosts = ["fisher", "alice2", "rex", "ficus", "batty", "localhost"]

parser = argparse.ArgumentParser(description='Slay sim clients on list of hosts.',
	epilog = 'slayClients.py [-a] [hosts], where no args is all default hosts, -a adds \
	default hosts to list, and hosts is the list of hosts, if no -a  and hosts, then use \
	only hosts')
parser.add_argument("-a", "--all", action='store_true', help="all hosts")
parser.add_argument("hosts", nargs=argparse.REMAINDER, help='list of hosts')
args = parser.parse_args()

if (args.all or args.hosts == []):
	hosts = allhosts + args.hosts
else:
	hosts = args.hosts

for h in hosts:
	output = check_output(['/Users/steve/bin/slayhosts.py', command, h])
	print(output.decode("utf-8"), end='')
