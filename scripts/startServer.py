#!/usr/bin/env python3

from subprocess import check_output
import argparse

job = "/Users/steve/sim/plasticity/02.sensitivity/scripts/jobRepeat.bash"
path = "/Users/steve/sim/grpcControl/sim_server"
control = "/Users/steve/sim/grpcControl/control/sensitivity.control"
host = "fisher"

parser = argparse.ArgumentParser(description='Start sim clients on list of hosts.')
parser.add_argument("host", nargs=argparse.REMAINDER, help='host for server')
args = parser.parse_args()

if not args.host == []:
	host = args.host[0]

check_output(['ssh', host, job, '1', path, control])
output = check_output(['ssh', host, 'cat', '/tmp/jobout*.out'])
print(output.decode("utf-8"), end='')
