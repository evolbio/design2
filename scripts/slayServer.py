#!/usr/bin/env python3

from subprocess import check_output
import argparse

command = "/Users/steve/sim/grpcControl/sim_shutdown"
server = "fisher"

parser = argparse.ArgumentParser(description='Shutdown server.',
	usage = 'slayServer.py [-h] -g | -n [server]',
	epilog = 'slayServer.py [-h] -g | -n [server], where no args returns help; \
	-g for graceful shutdown; -n for now shutdown; \
	optional server, default is fisher')
parser.add_argument("-g", "--graceful", action='store_true', help="graceful shutdown")
parser.add_argument("-n", "--now", action='store_true', help="now shutdown")
parser.add_argument("server", nargs=argparse.REMAINDER, help='server name')
args = parser.parse_args()

if (not args.now and not args.graceful) or (args.now and args.graceful):
	parser.print_help()
	print()
	exit()
	
if not args.server == []:
	server = args.server[0]
	
if args.now:
	now = 'now'
else:
	now = ''

output = check_output(['ssh', server, command, now])
print(output.decode("utf-8"), end='')

output = check_output(['/Users/steve/sim/plasticity/02.sensitivity/scripts/statusServer.py', server])
print(output.decode("utf-8"), end='')
