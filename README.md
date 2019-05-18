# GENIE-SMEM
MIT 6.890 Final Project - Using ML to predict suffix array location to discover SMEMs

Our implementations use the BWA-SMEM, LUT and RMI approaches to finding SMEM's within large query sequences. You may change the reference sequence and query sequence for each approach.

Requirements: Python 3, biopython (pip install biopython)

Exact Match
1) cd into GENIE/SMEM
2) run python ExactMatch.py


Super Maximal Exact Match (SMEM)
  - cd into GENIE/SMEM 
  - run python SMEM.py to run all three types of SMEM search: BWA, LUT, and RMI
    
   1) LUT
      - Run python LUT.py to create a new LUT with the given parameters. The LUT size can be modified at the end of the file. 
    
   2) RMI
      - Run python RMI_LUT.py to create a new RMI with the given parameters. The RMI size and expert levels can be modified at         the end of the file.
