# EICreconOutputReader

Contact - Kong Tu (kongtu@bnl.gov)

This is a branch for benchmark scripts to analyze EIC reconstructed data or MC. The input of this reader is the EICrecon Output, as the name suggested. I have two options: i) RDataFrame used previously in the physics benchmarks in ATHENA and ii) the TTree reader. However, for simplicity, here only shows option ii).

This repo will be actively developed and added with more examples as we progress, e.g., how to analyze a DIS physics event or an exclusive physics event. Hopefully, there will be more examples showing up under `src`. Contributions are welcomed. 


## Setting up the environment

Fresh start with an eic_shell, if one hasn't done it yet:

```
wget --output-document install.sh http://get.epic-eic.org --no-check-certificate
	
bash install.sh
```

Run it:

```./eic-shell```

## EICrecon reader

Install EICreconOutputReader

```git clone https://github.com/KongTu/EICreconOutputReader.git --branch benchmark-july-2023```

## Running diffractive VM (e.g., phi) samples:

Previously, we access the files via S3.
Look into getInputFromS3-exclusive.sh to modify accordingly what to grab from S3:

```./getInputFromS3-exclusive.sh``` (this is grabing a few files on s3 for an example)

Note that for official compaign simulation, the software simulation or validation team knows better where the generators are. Simply replace `input/rec-batch_5_official_\*.eicrecon.tree.edm4eic.root` to whatever the directory+name will be.

Now, it is much easier to use xrootd. So replace the 1st argument with 
```root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/23.12.0/epic_craterlake/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/Sartre/Coherent/sartre_bnonsat_Au_phi_ab_eAu_1.0000.eicrecon.tree.edm4eic.root```, where one can change the file name or directory. Wildcard * is supported.

Run with TTreeReader:

```./runDiffractiveVMReader.sh input/rec-batch_5_official_\*.eicrecon.tree.edm4eic.root output/eicrecon-sartre_coherent_phi```

with input files in `input`, output results in `output`, and benchmark figures in `figures`.
