#import dependencies
import os
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
class NUMT_finder():
	def __init__(self, organism_name):
		self.organism_name=str(organism_name)

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
				current_assembly=ftp_path.split('/')[-1]
				#download latest assembly genome
				call(
						f'wget --output-document=../data/gDNA.fna.gz {ftp_path}/{current_assembly}_genomic.fna.gz',
						shell=True
					)
				#decompress gDNA
				call(f'gzip -d ../data/gDNA.fna.gz',shell=True)
		except:
			print(f'A problem occured during {self.organism_name} gDNA acquisition!')

	def LASTALignment(self):
		try:
			if (os.path.getsize('../data/mtDNA.fna')>100):
				#generate LASTAL db
				call(
						'lastdb ../data/db ../data/gDNA.fna',
						shell=True
					)
				#align gDNA with mtDNA
				call(
						f"""lastal -r{lastal_settings['match_score']} -q{lastal_settings['mismatch_score']} -a{lastal_settings['gap_opening_score']} -b{lastal_settings['gap_extension_score']} ../data/db ../data/dmtDNA.fna  > ../data/aligned_dmtDNA.afa""",
						shell=True
					)
			else:
				print(f'A problem occured during {self.organism_name} db building or alignment!')
		except:
			print(f'A problem occured during {self.organism_name} db building or alignment!')

	def process_alignment(self):
		try:
			df_input=[]
			#get data from alignment file
			with open('../data//aligned_dmtDNA.afa')as infile:
				content=infile.readlines()
				for index, line in enumerate(content):
					if 'score' in line:
						#general information
						score,eg2_value,e_value=int(line.rsplit()[1].split('=')[1]),float(line.rsplit()[2].split('=')[1]),float(line.rsplit()[3].split('=')[1])
						#genomic information
						genomic=content[index + 1]
						genomic_id,genomic_start,genomic_length,genomic_strand,genomic_size,genomic_sequence=genomic.rsplit()[1],int(genomic.rsplit()[2]),int(genomic.rsplit()[3]),genomic.rsplit()[4],int(genomic.rsplit()[5]),genomic.rsplit()[6]
						#mitochondrial information
						mitochondrial=content[index + 2]
						mitochondrial_start,mitochondrial_length,mitochondrial_strand,mitochondrial_sequence = int(mitochondrial.rsplit()[2]),int(mitochondrial.rsplit()[3]),mitochondrial.rsplit()[4],mitochondrial.rsplit()[6]
						df_input.append([
								score, eg2_value, e_value, genomic_id, genomic_start,mitochondrial_start,
								genomic_length, mitochondrial_length, genomic_strand,mitochondrial_strand,
								genomic_size, genomic_sequence,mitochondrial_sequence
							])
			#create df from alignment info
			alignments=pd.DataFrame(
					data=df_input,
					columns=[
						'score', 'eg2_value', 'e_value', 'genomic_id', 'genomic_start',
						'mitochondrial_start', 'genomic_length', 'mitochondrial_length', 'genomic_strand',
						'mitochondrial_strand', 'genomic_size', 'genomic_sequence',
						'mitochondrial_sequence'
					]
				)
			#create filter to discard artifacts that are the results of using double mtDNA
			mtRecord=SeqIO.read("../data/mtDNA.fna", "fasta")
			mtSize=len(str(mtRecord.seq))
			size_fil=alignments['mitochondrial_start']<mtSize
			#get current assembly for filename
			#get latest assembly
			ftp_path=assembly_summary.loc[
				assembly_summary['organism_name']==self.organism_name.replace('_',' ').capitalize()
			]['ftp_path'].tolist()[0]
			current_assembly=ftp_path.split('/')[-1]
			#apply filters and write output into results folder
			alignments=alignments[size_fil]
			alignments.to_csv(f'../results/{self.organism_name}_{current_assembly}_numts.csv',header=True)
			#remove unnecessary file
			call("""rm ../data/* !('mitochondrion.1.1.genomic.fna'|'assembly_summary_refseq.txt'|'organism_names.txt'|'settings.json')""",shell=True)
		except:
			print(f'A problem occured during {self.organism_name} alignment processing!')

#download mitochondrial file
if os.path.exists('../data/mitochondrion.1.1.genomic.fna')==False:
	call(
			f'wget --directory-prefix=../data/ https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz',
			shell=True
		)
	call(
			f'gzip -d ../data/mitochondrion.1.1.genomic.fna.gz',
			shell=True
		)

if organism_names.size>1:
	for organism_name in organism_names:
		NUMT_class=NUMT_finder(organism_name)
		NUMT_class.get_mtDNA()
		NUMT_class.get_gDNA()
		NUMT_class.LASTALignment()
		NUMT_class.process_alignment()
else:
	NUMT_class=NUMT_finder(organism_names)
	NUMT_class.get_mtDNA()
	NUMT_class.get_gDNA()
	NUMT_class.LASTALignment()
	NUMT_class.process_alignment()