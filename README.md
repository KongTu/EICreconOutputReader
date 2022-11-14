# EICreconOutputReader

Contact - Kong Tu (kongtu@bnl.gov)

This is a reader script to analyze EIC reconstructed data or MC. The input of this reader is the EICrecon Output, as the name suggested. Here I took the RDataFrame used previously in the physics benchmarks. 

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

```git clone https://github.com/KongTu/EICreconOutputReader.git```

## Running samples:

Look into getInputFromS3.sh to modify accordingly what to grab from S3:

```./getInputFromS3.sh```

Run the singleParticleReader:

```./runSingleParticleReader.sh input/INPUT_NAME.root OUTPUT_NAME```
