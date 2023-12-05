# EICreconOutputReader

Contact - Kong Tu (kongtu@bnl.gov)

This is a reader script to analyze EIC reconstructed data or MC. The input of this reader is the EICrecon Output, as the name suggested. Here I took the RDataFrame used previously in the physics benchmarks. 

This repo will be actively developed and added with more examples as we progress, e.g., how to analyze a DIS physics event or an exclusive physics event. Hopefully, there will be more examples showing up under `src`, e.g., `src/readSingleParticles.cxx`. Contributions are welcomed. 

In addition, currently this code can utilize pfRICH-configs to study different designs/optimatizations of the backward pfRICH PID capbability.

## Setting up the environment

Fresh start with an eic_shell, if one hasn't done it yet:

```
wget --output-document install.sh http://get.epic-eic.org --no-check-certificate
	
bash install.sh
```

Run it:

```./eic-shell```

Set your LD_LIBRARY_PATH	

```export LD_LIBRARY_PATH=/gpfs02/eic/ztu/EPIC/software_tutorial/analysis/lib:${LD_LIBRARY_PATH}```

## Additional install/setup if one doesn't want to point to Kong's tmp directory...

For this example, I am going to set a soft link:

```cd /tmp```

```ln -s <your-directory-where-the-eic-shell-is> EPIC-Kong```

Install the data model - EDM4EIC with irt data model branch

```
git clone https://github.com/eic/EDM4eic.git --branch irt-data-model
cd EDM4eic && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/tmp/EPIC-Kong ..
make -j8 install
```

Install the IRT

```
git clone https://github.com/eic/irt.git
cd irt && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/tmp/EPIC-Kong ..
make -j8 install
```

## EICrecon reader

Install EICreconOutputReader

```git clone https://github.com/KongTu/EICreconOutputReader.git```

Note that under `src/pleaseIncludeMe.h`, there is a direct include: 

`#include "/gpfs02/eic/ztu/EPIC/software_tutorial/analysis/irt/delphes/include/DelphesConfig.h"`. 

For the time being, one can just leave it as it is.

## Running samples:

Look into getInputFromS3.sh to modify accordingly what to grab from S3:

```./getInputFromS3.sh```

For example, for DIS events from the Oct Simulation Campaign. 

```./getInputFromS3-DIS.sh 18 275```
the 1st and 2nd arguments are energy configuration of DIS.

Run the singleParticleReader:

```./runSingleParticleReader.sh input/INPUT_NAME.root OUTPUT_NAME```


