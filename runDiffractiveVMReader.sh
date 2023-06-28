#!/bin/bash

input=$1
output=$2

echo "Running analysis"
echo "Input edm4eic data is at [ ${input} ] "
echo "Output at [ ${output} ]"

root -b -q src/diffractive_vm_simple_analysis.cxx+\(\"${input}\",\"${output}\"\)

# plotting benchmark figures.
# figures are stored under figures.
root -b -q macros/plot_diffractive_vm_physics_benchmark.C\(\"${output}\"\)
root -b -q macros/plot_diffractive_vm_resolution.C\(\"${output}\"\)
root -b -q macros/plot_diffractive_event_kinematics.C\(\"${output}\"\)
