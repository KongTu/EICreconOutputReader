#!/bin/bash

mc config host add S3 https://dtn01.sdcc.bnl.gov:9000 eicS3read eicS3read

for i in {1..9}
do
 	mc cp --insecure S3/eictest/EPIC/RECO/22.11.2/epic_brycecanyon/DIS/NC/10x100/minQ2=1/pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_vtxfix_1.000${i}.eicrecon.tree.edm4eic.root ./input/input-eicrecon-DIS_10x100_Q2-1_000${i}.root

done

for i in {10..99}
do
    mc cp --insecure S3/eictest/EPIC/RECO/22.11.2/epic_brycecanyon/DIS/NC/10x100/minQ2=1/pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_vtxfix_1.00${i}.eicrecon.tree.edm4eic.root ./input/input-eicrecon-DIS_10x100_Q2-1_00${i}.root

done

for i in {100..999}
do
    mc cp --insecure S3/eictest/EPIC/RECO/22.11.2/epic_brycecanyon/DIS/NC/10x100/minQ2=1/pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_vtxfix_1.0${i}.eicrecon.tree.edm4eic.root ./input/input-eicrecon-DIS_10x100_Q2-1_0${i}.root

done

echo "************************************************"
echo "Download from S3 successful."
echo "Next step is to run reader on them, for example:"
echo "************************************************"
echo "./runSingleParticleReader.sh input/input-eicrecon-DIS_10x100_Q2-1_\*.root eicrecon-DIS_10x100_Q2-1"
echo "************************************************"

