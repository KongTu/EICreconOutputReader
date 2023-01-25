#!/bin/bash

input=$1
output=$2
VM=$3

echo "Running analysis"
echo "Input edm4eic data is at [ ${input} ] "
echo "Output at [ ${output} ]"
echo "VM species is [ ${VM} ]"

root -b -q src/diffractive_vm_simple_analysis.cxx+\(\"${input}\",\"${output}\",${VM}\)


