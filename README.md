# NUMT_finder
 Application for finding NUMTs in NCBI reference genomes.

 # Upon usage, please cite
 ´´´bibtex
 @article{biro2024mitochondrial,
  title={Mitochondrial genome plasticity of mammalian species},
  author={Bir{\'o}, B{\'a}lint and G{\'a}l, Zolt{\'a}n and Fekete, Zs{\'o}fia and Klecska, Eszter and Hoffmann, Orsolya Ivett},
  journal={BMC Genomics},
  volume={25},
  number={1},
  pages={1--14},
  year={2024},
  publisher={BioMed Central}
}
 ´´´

 # Abstract
There is an ongoing process in which mitochondrial sequences are being integrated into the nuclear genome. The importance of these sequences has already been revealed in cancer biology, forensic, phylogenetic studies and in the evolution of the eukaryotic genetic information. Human and numerous model organisms’ genomes were described from those sequences point of view. Furthermore, recent studies were published on the patterns of these nuclear localised mitochondrial sequences in different taxa.

However, the results of the previously released studies are difficult to compare due to the lack of standardised methods and/or using few numbers of genomes. Therefore, in this paper our primary goal is to establish a uniform mining pipeline to explore these nuclear localised mitochondrial sequences.

Our results show that the frequency of several repetitive elements is higher in the flanking regions of these sequences than expected. A machine learning model reveals that the flanking regions' repetitive elements and different structural characteristics are highly influential during the integration process.

In this paper, we introduce a general mining pipeline for all mammalian genomes. The workflow is publicly available and is believed to serve as a validated baseline for future research in this field. We confirm the widespread opinion, on - as to our current knowledge - the largest dataset, that structural circumstances and events corresponding to repetitive elements are highly significant. An accurate model has also been trained to predict these sequences and their corresponding flanking regions

![graphical_abstract](/data/fig6.png)

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
```txt
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

Used Python packages (if the module is not built in, the version number and conda installation are provided):
---
- os
- json
- ftplib
- subprocess
- numpy (1.21.6) `conda install numpy`
- pandas (1.3.5) `conda install -c anaconda pandas`
- Bio (1.78) `conda install -c conda-forge biopython`
