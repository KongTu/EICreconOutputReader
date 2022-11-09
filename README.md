# EICreconOutputReader

contact - Kong Tu (kongtu@bnl.gov)

This is a reader code to analyze EIC reconstructed data or MC. The detail instructions will be documented along the development. 

One can use pfRICH-configs to study backward PID capbability.

To get started:
- look into getInputFromS3.sh to modify accordingly what to grab from S3;
```./getInputFromS3.sh```

- Run the singleParticleReader:
```./runSingleParticleReader.sh input/INPUT_NAME.root OUTPUT_NAME```
