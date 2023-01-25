# NUMT_finder
 Application for finding NUMTs in NCBI reference genomes.

 At this moment, the program is applicable in Unix based systems. All the codes are written in Python (3.7.13). By default, the data directory contains seetings.json with the corresponding LASTAL settings; and organism_names.txt file with the organism names.

 Setting up the environment:
---
```bash
mkdir code data results
conda install -c bioconda last 
conda install -c bioconda samtools (1.6)
```

Used Python packages and external programs (if the module is not built in, the version number and conda installation are provided):
---
- os
- json
- ftplib
- subprocess
- numpy (1.21.6) `conda install numpy`
- pandas (1.3.5) `conda install -c anaconda pandas`
- Bio (1.78) `conda install -c conda-forge biopython`