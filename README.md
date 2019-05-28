# GENIE-SMEM
MIT 6.890 Final Project - Using ML to predict suffix array location to discover SMEMs

Our implementations use the BWA-SMEM, LUT and RMI approaches to finding SMEM's within large query sequences, as described in the paper
in the docs directory.

Requirements: Python 3, biopython (pip install biopython), as well as numpy and scipy


## Overview

An overview of the workflow, where things are kept, and the data structures in this repo. 

### Data

The data is kept in the [SMEM/data](SMEM/data) directory. This consists of the actual reference 
sequence, the generated FM indexes, the generated LUTs, and any queries used for exact match.

### Backwards Search Algorithm

The backwards search algorithm is implemented in [SMEM/ExactMatch.py](SMEM/ExactMatch.py) in the exact_match_back_prop 
method which returns suffix array location or the exact_match method which returns the actual 
postion in the reference sequence. 

### Exact Match

This class, [SMEM/ExactMatch.py](SMEM/ExactMatch.py),  is able to create FM indexes, store and load them, do exact match via the backwards 
search algorithm, and generate and save queries from pieces of the reference sequence. Importantly, this uses an unoptimized version 
of the BW algorithm, and so when creating the BW Matrix it requires n^2 space where n is the size of the reference sequence. 
This means that for this code it is impractical to create an FM index for a ref sequence of much more than 
100,000 bases. To use a larger data set, use an optimized version of the BW algorithm to create the FM index.

### Look up Table

This class, [SMEM/LUT.py](SMEM/LUT.py), is used to create and store a look up table given 
a key size K, and an ExactMatch instance, (used to get the ref sequence and suffix array). The lookup 
table maps every subsequence of size K to its suffix array position. 

### RMI 

The actual code used to train and predict for the RMI is located in [SMEM/RMI.py](SMEM/RMI.py). We built a class 
around this in [SMEM/RMI_LUT.py](SMEM/RMI_LUT.py) which given an expert level, a data file, and a key size K, trains 
an RMI on the sequence in the data file given which takes as input a string of size K and returns the suffix position of that 
string. It will return a higher lower bound then upper bound if the sequence is not contained in 
the reference sequence. It is able to be pickled and saved as a .pkl file. 

### SMEM

We have implemented three ways to find SMEMs in the [SMEM/SMEM.py](SMEM/SMEM.py) class, using the traditional 
BWA-SMEM algorithm, using an LUT-SMEM algorithm, and using an RMI-SMEM algorithm. The workings of these algorithms are
described in the writeup in the [docs](docs/) directory. Simply pass the SMEM class 
an instance of the ExactMatch class (after loading the fm index) and then choose which 
type of SMEM search you want to run.



## Running 

### Exact match

To run the exact match algorithm with backwards search, use the ExactMatch.py class. To run the exact match algorithm with the RMI
run the train.py script. 


### SMEM

To run any of the three SMEM algorithms, use the SMEM.py class. Pass in an instance of the ExactMatch and 
then choose which type of SMEM algorithm you want, and pass it a query. There is a method for 
creating queries randomly and from creating them from sections of the reference sequence. 
For optimal performance, the key size for the LUT or RMI must be close to that of the predicted 
average SMEM size. 


## Info

For more information, first read the paper in docs, and then feel free to contact either of the authors (contact info in the paper)!

For a different implementation of exact match, see https://github.com/ashwatht/GENIE, which has a combination of c++, c, and 
python implementations. Much of the work in this repo is adapted from the ideas or code used originally in the GENIE project and paper.
