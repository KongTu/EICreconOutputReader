# EICreconOutputReader

Contact - Kong Tu (kongtu@bnl.gov)

This is a reader script to analyze EIC reconstructed data or MC. The input of this reader is the EICrecon Output, as the name suggested. Here I took the i) RDataFrame used previously in the physics benchmarks in ATHENA and ii) the TTree reader. 

This repo will be actively developed and added with more examples as we progress, e.g., how to analyze a DIS physics event or an exclusive physics event. Hopefully, there will be more examples showing up under `src`, e.g., `src/readSingleParticles.cxx`. Contributions are welcomed. 


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

```git clone https://github.com/KongTu/EICreconOutputReader.git --branch simple-branch```

## Running single particle samples:

Look into getInputFromS3.sh to modify accordingly what to grab from S3:

```./getInputFromS3.sh```

Run with RDataFrame:

```./runSingleParticleReader.sh input/INPUT_NAME.root OUTPUT_NAME```

Or Run with TTreeReader:

```./runSingleParticleSimpleReader.sh input/INPUT_NAME.root OUTPUT_NAME```

## Running diffractive VM (e.g., phi) samples:

Look into getInputFromS3-exclusive.sh to modify accordingly what to grab from S3:

```./getInputFromS3-exclusive.sh``` (its grabing 1000 files on s3 in this example, give it a few mins)

Run with TTreeReader:

```./runDiffractiveVMReader.sh input/rec-batch_5_official_\*.eicrecon.tree.edm4eic.root output/eicrecon-sartre_coherent_phi```

with input files in `input` and output results in `output`
