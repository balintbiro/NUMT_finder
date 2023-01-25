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
assembly_summary=pd.read_csv('../data/assembly_summary_refseq.txt',sep='\t',skiprows=1)

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
		self.organism_name=str(organism_name)
		self.current_assembly=''

	def get_mtDNA(self):
		try:
			#get mitochondrial id
			mtID=run(
					f"""egrep '{self.organism_name.replace('_',' ').capitalize()}' ../data/mitochondrion.1.1.genomic.fna | grep mitochondrion""",
					shell=True,capture_output=True
				)
			mtID=str(mtID.stdout).split()[0][3:]
			#get mitochondrial sequence
			call(
					f'samtools faidx ../data/mitochondrion.1.1.genomic.fna {mtID} > ../data/mtDNA.fna',
					shell=True
				)
			#get duplicated mitochondria
			mtRecord=SeqIO.read("../data/mtDNA.fna", "fasta")
			mtSeq=str(mtRecord.seq)
			with open('../data/dmtDNA.fna','w')as outfile:
				outfile.write('>'+str(mtID)+'\n'+2*mtSeq)
		except:
			print(f'A problem occured during {self.organism_name} mtDNA acquisition!')

	def get_gDNA(self):
		try:
			if os.path.exists('../data/mtDNA.fna')==False:
				print(f'No mtDNA sequence was found for {self.organism_name}!')
			elif os.path.exists('../data/gDNA.fna')==False:
				#get latest assembly
				ftp_path=assembly_summary.loc[
					assembly_summary['organism_name']==self.organism_name.replace('_',' ').capitalize()
				]['ftp_path'].tolist()[0]
				self.current_assembly=ftp_path.split('/')[-1]
				#download latest assembly genome
				call(
						f'wget --output-document=../data/gDNA.fna.gz {ftp_path}/{self.current_assembly}_genomic.fna.gz',
						shell=True
					)
				#decompress gDNA
				call(f'gzip -d ../data/gDNA.fna.gz',shell=True)
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
						../data/db ../data/dmtDNA.fna  > ../data/aligned_dmtDNA.afa""",
						shell=True
					)
				#remove unnecessary file
				call("""rm !('../data/mitochondrion.1.1.genomic.fna'|'../data/assembly_summary_refseq.txt|../data/organism_names.txt'|'../data/settings.json')""",shell=True)
			else:
				print(f'A problem occured during {self.organism_name} db building or alignment!')
		except:
			print(f'A problem occured during {self.organism_name} db building or alignment!')

#download mitochondrial file
if os.path.exists('../data/mitochondrion.1.1.genomic.fna')==False:#uncompressed file
	call(
			f'wget --directory-prefix=../data/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz',
			shell=True
		)
	call(
			f'gzip -d ../data/{filename}',
			shell=True
		)

if organism_names.size>1:
	for organism_name in organism_names:
		NUMT_class=get_genomes(organism_name)
		NUMT_class.get_mtDNA()
		NUMT_class.get_gDNA()
		NUMT_class.LASTALignment()
else:
	NUMT_class=get_genomes(organism_names)
	NUMT_class.get_mtDNA()
	NUMT_class.get_gDNA()
	NUMT_class.LASTALignment()