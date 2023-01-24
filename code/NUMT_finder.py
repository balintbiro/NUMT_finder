#import dependencies
import os
import re
import json
import numpy as np
import pandas as pd
from Bio import SeqIO
from ftplib import FTP
from subprocess import call, run

#get organism names
organism_names=np.loadtxt('../data/organism_names.txt',dtype=str,delimiter='\n')

#get assembly summary file
if os.path.exists('../data/assembly_summary_refseq.txt')==False:
	ftp_site='https://ftp.ncbi.nlm.nih.gov/genomes/refseq/'
	filename='assembly_summary_refseq.txt'
	call(f'wget --directory-prefix=../data/ {ftp_site}{filename}',shell=True)

#read assembly summary file
assembly_summary=pd.read_csv(
		'../data/assembly_summary_refseq.txt',
		sep='\t',skiprows=1
	)

#sort assembly summary df based on genome quality
assembly_qualities=['Complete Genome','Chromosome','Scaffold','Contig']
assembly_summary['assembly_quality']=pd.Categorical(
		assembly_summary['assembly_level'],
		categories=assembly_qualities,ordered=True
	)
assembly_summary=assembly_summary.sort_values(by=['assembly_quality'])

#read LASTAL setting .json file
with open('../data/settings.json')as infile:
	lastal_settings=infile.read()
lastal_settings=json.loads(lastal_settings)

#class definition for processing DNA
class get_genomes():
	def __init__(self, organism_name):
		self.organism_name=organism_name

	def get_mtDNA(self):
		try:
			#get mitochondrial id
			mtID=run(
					f"""egrep '{self.organism_name.replace('_',' ')}' ../data/genomes/mitochondrion.1.1.genomic.fna | grep mitochondrion""",
					shell=True,capture_output=True
				)
			mtID=str(mtID.stdout).split()[0][3:]
			#get mitochondrial sequence
			call(
					f'samtools faidx ../data/genomes/mitochondrion.1.1.genomic.fna {mtID} > ../data/genomes/mtDNA.fna',
					shell=True
				)
			#some organisms are not present in the mitochondrion.1.1.genomic.fna file so try mitochondrion.2.1.genomic.fna instead
			if os.path.getsize('../data/genomes/mtDNA.fna')<1000:
				mtID=run(
						f"""egrep '{self.organism_name.replace('_',' ')}' ../data/genomes/mitochondrion.2.1.genomic.fna | grep mitochondrion""",
						shell=True,capture_output=True
					)
				mtID=str(mtID.stdout).split()[0][3:]
				call(
						f'samtools faidx ../data/genomes/mitochondrion.2.1.genomic.fna {mtID} > ../data/genomes/mtDNA.fna',
						shell=True
					)
			#get duplicated mitochondria
			mtRecord=SeqIO.read("../data/genomes/mtDNA.fna", "fasta")
			mtSeq=str(mtRecord.seq)
			with open('../data/genomes/dmtDNA.fna','w')as outfile:
				outfile.write('>'+str(mtID)+'\n'+2*mtSeq)
		except:
			print(f'A problem occured during {self.organism_name} mtDNA acquisition!')

	def get_gDNA(self):
		try:
			if os.path.exists('../data/mtDNA.fna')==False:
				print(f'No mtDNA sequence was found for {self.organism_name}!')
			else:
				#get latest assembly
				ftp_path=assembly_summary.loc[assembly_summary['organism_name'].str.contains(organism_name.replace('_',' '),flags=re.IGNORECASE)]['ftp_path'].tolist()[0]
				current_assembly=ftp_path.split('/')[-1]
				#download latest assembly genome
				call(
						f'wget --output-document=../data/gDNA.fna.gz {ftp_path}/{current_assembly}_genomic.fna.gz',
						shell=True
					)
				#decompress gDNA
				call(f'gzip -d ../data/gDNA.fna.gz',shell=True)
				print(f'For {self.organism_name} the current assembly ID is: {current_assembly}!')
		except:
			print(f'A problem occured during {self.organism_name} gDNA acquisition!')

	def LASTALignment(self):
		try:
			if os.path.getsize('../data/mtDNA.fna')>100:
				#generate LASTAL db
				call(
						'lastdb ../data/db ../data/gDNA.fna',
						shell=True
					)
				#align gDNA with mtDNA
				call(
						f"""lastal -r{lastal_settings['match_score']}
						-q{lastal_settings['mismatch_score']}
						-a{lastal_settings['gap_opening_score']}
						-b{lastal_settings['gap_extension_score']}
						../data/genomes/db ../data/genomes/dmtDNA.fna  > ../data/genomes/aligned_dmtDNA.afa""",
						shell=True
					)
			else:
				print(f'A problem occured during {self.organism_name} db building or alignment!')
		except:
			print(f'A problem occured during {self.organism_name} db building or alignment!')

#function to download mitochondrial files
def get_all_mtDNA(filename):
	call(f'wget --directory-prefix=../data/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/{filename}')
	call(f'gzip -d ../data/{filename}')

#download mitochondrial files
mitochondrial_files=pd.Series(['mitochondrion.1.1.genomic.fna.gz','mitochondrion.2.1.genomic.fna.gz'])
mitochondrial_files.apply(get_all_mtDNA)

for organism_name in organism_names:
	NUMT_class=get_genomes(organism_name)
	NUMT_class.get_mtDNA()
	NUMT_class.get_gDNA()
	NUMT_class.LASTALignment()