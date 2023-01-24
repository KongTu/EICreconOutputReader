mc config host add S3 https://dtn01.sdcc.bnl.gov:9000 eicS3read eicS3read

echo "Download simulation files from S3 for exclusive ${1},${2} eA coherent events for 18x110 GeV configuration Q2>1"

for i in {1..9}
do
 	mc cp --insecure S3/eictest/EPIC/RECO/22.11.3/epic_arches/EXCLUSIVE/DIFFRACTIVE_${1}_ABCONV/Sartre/Coherent/sartre_bnonsat_Au_${2}_ab_eAu_1.000${i}.eicrecon.tree.edm4eic.root ./input/rec-batch_5_${2}_official_000${i}.eicrecon.tree.edm4eic.root
done

for i in {10..99}
do
    mc cp --insecure S3/eictest/EPIC/RECO/22.11.3/epic_arches/EXCLUSIVE/DIFFRACTIVE_${1}_ABCONV/Sartre/Coherent/sartre_bnonsat_Au_${2}_ab_eAu_1.00${i}.eicrecon.tree.edm4eic.root ./input/rec-batch_5_${2}_official_00${i}.eicrecon.tree.edm4eic.root
done

for i in {100..999}
do
    mc cp --insecure S3/eictest/EPIC/RECO/22.11.3/epic_arches/EXCLUSIVE/DIFFRACTIVE_${1}_ABCONV/Sartre/Coherent/sartre_bnonsat_Au_${2}_ab_eAu_1.0${i}.eicrecon.tree.edm4eic.root ./input/rec-batch_5_${2}_official_0${i}.eicrecon.tree.edm4eic.root

done
