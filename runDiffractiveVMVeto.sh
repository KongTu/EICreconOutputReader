#!/bin/bash

input=$1
output=$2
veto=$3

echo "Running analysis"
echo "Input edm4eic data is at [ ${input} ] "
echo "Output at [ ${output} ]"
echo "Are we vetoing? " ${veto}

root -b -q src/diffractive_vm_veto_analysis.cxx+\(\"${input}\",\"${output}\"\,${veto}\)