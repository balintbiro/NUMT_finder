# NUMT_finder
 Application for finding NUMTs in NCBI reference genomes.

 At this moment, the program is applicable in Unix based systems. All the codes are written in Python (3.7.13). By default, the data directory contains seetings.json with the corresponding LASTAL settings; and organism_names.txt file with the organism names.

The content of settings.json file (default settings of lastal:
---
```json
{
	"match_score":1,
	"mismatch_score":1,
	"gap_opening_score":7,
	"gap_extension_score":1
}
```

The organism_names.txt file contains the organism names line by line, all lowercase without comma as the following:
```bash
mus musculus
homo sapiens
```

 Setting up the environment:
---
```bash
mkdir code data results
conda install -c bioconda last 
conda install -c bioconda samtools (1.6)
```

External program usage example.
---
The external programs (lasta and samtools) and command line tools (grep, egrep, wget, rm etc.) were run from Python with the subprocess module as followings:
```python
from subprocess import call, run

call(f'gzip -d ../data/gDNA.fna.gz',shell=True)
mtID=run(
		f"""egrep '{self.organism_name.replace('_',' ').capitalize()}'
		../data/mitochondrion.1.1.genomic.fna | grep mitochondrion""",
		shell=True,capture_output=True
		)
```

Used Python packages (if the module is not built in, the version number and conda installation are provided):
---
- os
- json
- ftplib
- subprocess
- numpy (1.21.6) `conda install numpy`
- pandas (1.3.5) `conda install -c anaconda pandas`
- Bio (1.78) `conda install -c conda-forge biopython`