#!/bin/bash

mc config host add S3 https://dtn01.sdcc.bnl.gov:9000 eicS3read eicS3read

echo "Download simulation files from S3 for DIS NC events for ${1}x${2} GeV configuration Q2>1"

mc cp --insecure S3/eictest/EPIC/RECO/22.11.2/epic_brycecanyon/DIS/NC/${1}x${2}/minQ2=1/pythia8NCDIS_${1}x${2}_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_vtxfix_1.0001.eicrecon.tree.edm4eic.root ./input/input-eicrecon-DIS_${1}x${2}_Q2-1_0001.root




echo "************************************************"
echo "Download from S3 successful."
echo "Next step is to run reader on them, for example:"
echo "************************************************"
echo "./runSingleParticleReader.sh input/input-eicrecon-DIS_${1}x${2}_Q2-1_\*.root eicrecon-DIS_${1}x${2}_Q2-1"
echo "************************************************"

