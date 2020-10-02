#!/usr/bin/env bash

linefile=../data/line_data

display_usage() {
  echo -e "\nUsage:\n$(basename $0) [inputFolder] [Subgraph] [LF] [VF] \n"
  echo -e "   where: Subgraph, LF and VF can either be 1 or 0\n"
}

if [ "$#" -ne 4 ]; then
  echo "Illegal number of parameters. Expected 4 parameters"
  display_usage
  exit 1
fi

for i in $1/*.dat; do
  ../bin/alg -u Cmdenv --cmdenv-express-mode=true -f <(echo "
[General]
network=Custom1
Custom1.node[**].line_filename=\"$linefile\"
Custom1.node[**].bus_filename=\"$i\"
Custom1.node[**].dodiffCommunication=$2
Custom1.node[**].zeroloadslack=$3
Custom1.node[**].zerovoltslack=$4
") ;done
