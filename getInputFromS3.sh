#!/bin/bash

mc config host add S3 https://dtn01.sdcc.bnl.gov:9000 eicS3read eicS3read
mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/pi-/1GeV/130to177deg/pi-_1GeV_130to177deg.0006.eicrecon.tree.edm4eic.root ./input/input-pi-1GeV-backward.root
mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/pi-/1GeV/3to50deg/pi-_1GeV_3to50deg.0006.eicrecon.tree.edm4eic.root ./input/input-pi-1GeV-forward.root
mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/pi-/1GeV/45to135deg/pi-_1GeV_45to135deg.0001.eicrecon.tree.edm4eic.root ./input/input-pi-1GeV-mid.root

mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/pi-/10GeV/130to177deg/pi-_10GeV_130to177deg.0006.eicrecon.tree.edm4eic.root ./input/input-pi-10GeV-backward.root
mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/pi-/10GeV/3to50deg/pi-_10GeV_3to50deg.0006.eicrecon.tree.edm4eic.root ./input/input-pi-10GeV-forward.root
mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/pi-/10GeV/45to135deg/pi-_10GeV_45to135deg.0001.eicrecon.tree.edm4eic.root ./input/input-pi-10GeV-mid.root

mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/e-/10GeV/130to177deg/e-_10GeV_130to177deg.0001.eicrecon.tree.edm4eic.root ./input/input-e-10GeV-backward.root
mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/e-/10GeV/3to50deg/e-_10GeV_3to50deg.0001.eicrecon.tree.edm4eic.root ./input/input-e-10GeV-forward.root
mc cp --insecure S3/eictest/EPIC/RECO/22.10.0/epic_brycecanyon/SINGLE/e-/10GeV/45to135deg/e-_10GeV_45to135deg.0001.eicrecon.tree.edm4eic.root ./input/input-e-10GeV-mid.root

echo "************************************************"
echo "Download from S3 successful."
echo "Next step is to run reader on them, for example:"
echo "************************************************"
echo "./runSingleParticleReader.sh input/input-pi-10GeV-backward.root pi-10GeV-backward"
echo "************************************************"
