# GENIE

## Initialization

On mac, start by running the following brew commands:
`brew install gnu-sed`
`brew install libomp`

## Setup 

Add a ```[DATA].fa``` file to the ```data/``` directory. For example, we have ```test.fa```. 

### Running on Linux

Replace the ```gsed``` commands in ```./genie.sh``` with ```sed```. You may also have to add an ```#include <string.h>``` to ```build_bwt_index/file.cpp```. 

## How to Run 
- Original implementations 
	- Naive: ```./genie.sh run_naive [DATA]```. For example, with ```test.fa```, run ```./genie.sh run_naive test```
	- Optimized: ```./genie.sh run_opt [DATA]```. For example, with ```test.fa```, run ```./genie.sh run_opt test```
- GENIE
	- Naive look up table (LUT): ```./genie.sh run_lut [DATA]```. For example, with ```test.fa```, run ```./genie.sh run_lut test```
	- RMI: ```./genie.sh run_rmi [DATA]```. For example, with ```test.fa```, run ```./genie.sh run_rmi test```

## Clean
Run ```./genie.sh clean [DATA]```. For example, with ```test.fa```, run ```./genie.sh clean test```.

