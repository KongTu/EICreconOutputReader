mc config host add S3 https://dtn01.sdcc.bnl.gov:9000 eicS3read eicS3read

echo "Download simulation files from S3 for exclusive phi eA coherent events for 18x110 GeV configuration Q2>1"

for i in {1..9}
do
 	mc cp --insecure S3/eictest/EPIC/RECO/22.11.3/epic_arches/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/Sartre/Coherent/sartre_bnonsat_Au_phi_ab_eAu_1.000${i}.eicrecon.tree.edm4eic.root ./input/rec-batch_5_official_000${i}.eicrecon.tree.edm4eic.root
done

# for i in {10..99}
# do
#     mc cp --insecure S3/eictest/EPIC/RECO/22.11.3/epic_arches/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/Sartre/Coherent/sartre_bnonsat_Au_phi_ab_eAu_1.00${i}.eicrecon.tree.edm4eic.root ./input/rec-batch_5_official_00${i}.eicrecon.tree.edm4eic.root
# done

# for i in {100..999}
# do
#     mc cp --insecure S3/eictest/EPIC/RECO/22.11.3/epic_arches/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/Sartre/Coherent/sartre_bnonsat_Au_phi_ab_eAu_1.0${i}.eicrecon.tree.edm4eic.root ./input/rec-batch_5_official_0${i}.eicrecon.tree.edm4eic.root
# done

echo "************************************************"
echo "Download from S3 successful."
echo "Next step is to run reader on them, for example:"
echo "************************************************"
echo "./runDiffractiveVMReader.sh input/rec-batch_5_official_\*.eicrecon.tree.edm4eic.root output/eicrecon-sartre_coherent_phi"
echo "************************************************"
