# EICreconOutputReader

contact - Kong Tu (kongtu@bnl.gov)

This is a reader code to analyze EIC reconstructed data or MC. The detail instructions will be documented along the development. 

One can use pfRICH-configs to study different designs/optimatizations of the backward pfRICH PID capbability.

## Setting up the environment

Fresh start with installing eic_shell:

```
	wget --output-document install.sh http://get.epic-eic.org --no-check-certificate
	
	bash install.sh
```

Run it:

```./eic-shell```

Set your LD_LIBRARY_PATH	

```export LD_LIBRARY_PATH=/tmp/EPIC-Kong/lib:${LD_LIBRARY_PATH}```

## Additional install/setup if one doesn't want to point to Kong's tmp dir...

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

`#include "/tmp/EPIC-Kong/irt/delphes/include/DelphesConfig.h"`. 

For the time being, one can just leave it as it is.

## Running samples:

Look into getInputFromS3.sh to modify accordingly what to grab from S3:

```./getInputFromS3.sh```

Run the singleParticleReader:

```./runSingleParticleReader.sh input/INPUT_NAME.root OUTPUT_NAME```
