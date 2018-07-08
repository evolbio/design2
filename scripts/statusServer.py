#!/usr/bin/env python3

from subprocess import check_output
import argparse

proc = 'sim_server'
command = ['ps', 'x', '-o', 'pid,%cpu,%mem,command', '|', 'grep']
host = ["fisher"]

parser = argparse.ArgumentParser(description='Status of server.')
parser.add_argument("host", nargs=argparse.REMAINDER, help='host, default is fisher')
args = parser.parse_args()

if not args.host == []:
	host = args.host
output = check_output(['ssh'] + host + command + [proc])
print("\n" + host[0] + ':')
for line in output.splitlines():
	if not b'statusServer' in line and not b'grep' in line:
		print("  " + line.decode("utf-8"))
print()