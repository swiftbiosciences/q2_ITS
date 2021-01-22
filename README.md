Swift’s ITS analysis workflow using Qiime 2 

Software package requirements:
Qiime 2 (https://qiime2.org/)

Setup:
1. Unzip the package file to get following files: q2wkflow_ITS_v3.sh,  
   mergeStrand_taxonomy.py, config.txt, primer file, and reference 
   file (1) for the 'sklearn' classifier using UNITE ITS reference 
   version 8. Please make sure to place file "mergeStrand_taxonomy.py" 
   and "q2wkflow_ITS_v3.sh" in the same file folder.

2. Edit “config.txt” to enter correct absolute paths to each tool, primer 
   file, expected read length after primer is trimmed, and clustering 
   similarity cutoff.

   Command to run:
   q2wkflow_ITS_v3.sh config.txt inputdir workdir

3. Results organization: ASV feature table can be found in "Export" directory; 
   "QObj" directory contains intermediate files and tables some of which could 
   be imported into Qiime2 web interface for further exploration and 
   visualization; "Fastq" directory contains intermediate fastq files at all
   different stages.


Additional notes:

1. This workflow is intended to serve as a basic example using Qiime 2
   for processing fungal ITS data prepared using Swift SNAP ITS1
   library kit. Users should feel free to make any modifications for
   improvements.

2. Additional quality trimming/filtering of reads may be needed depending on
   specific sequencing data sets.

3. The performance of this workflow is highly dependent on the taxonomy
   classifier chosen (this example uses naive baysian classifier
   with UNITE ITS reference set version 8 included in the package, for ITS
   analysis). Users should feel free to obtain and test  other type  of
   classifiers & reference sets that best fit the environments sampled.
