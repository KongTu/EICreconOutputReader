#!/bin/bash

mc config host add S3 https://dtn01.sdcc.bnl.gov:9000 eicS3read eicS3read
mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/pi-/10GeV/130to177deg/pi-_10GeV_130to177deg.0006.eicrecon.tree.edm4eic.root ./input/input-eicrecon-pi-10GeV-backward.root
mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/pi-/10GeV/130to177deg/pi-_10GeV_130to177deg.0006.juggler.tree.edm4eic.root ./input/input-juggler-pi-10GeV-backward.root

echo "************************************************"
echo "Download from S3 successful."
echo "Next step is to run reader on them, for example:"
echo "************************************************"
echo "./runSingleParticleReader.sh input/input-pi-10GeV-backward.root pi-10GeV-backward"
echo "************************************************"

