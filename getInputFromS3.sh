#!/bin/bash
mc config host add S3 https://dtn01.sdcc.bnl.gov:9000 eicS3read eicS3read
mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_arches/SINGLE/pi-/1GeV/130to177deg/pi-_1GeV_130to177deg.0006.eicrecon.tree.edm4eic.root ./input/input.root

echo "Next step is to run:"
echo "./runSingleParticleReader.sh input/input.root test"
