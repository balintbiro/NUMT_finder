#import dependencies
import os
import time
import json
import logging
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from ftplib import FTP
from subprocess import call, run
from argparse import RawTextHelpFormatter

class NUMT_finder:
	"""
	Class for getting NUMTs from mammalian reference genomes.
	"""
	def __init__(self, organism_name:str,out_path:str,settings:dict,assembly_summary:pd.DataFrame)->None:
		self.organism_name=str(organism_name)
		self.out_path=str(out_path)
		self.settings=settings
		self.assembly_summary=assembly_summary

	def get_mtDNA(self)->None:
		try:
			#get mitochondrial id
			mtID=run(
					f"""egrep '{self.organism_name.replace('_',' ').capitalize()}' {self.out_path}mitochondrion.1.1.genomic.fna | grep mitochondrion""",
					shell=True,capture_output=True
				)
			mtID=str(mtID.stdout).split()[0][3:]
			#get mitochondrial sequence
			call(
					f'samtools faidx {self.out_path}mitochondrion.1.1.genomic.fna {mtID} > {self.out_path}mtDNA.fna',
					shell=True
				)
			#get duplicated mitochondria
			mtRecord=SeqIO.read(f"{self.out_path}mtDNA.fna", "fasta")
			mtSeq=str(mtRecord.seq)
			with open(f'{self.out_path}dmtDNA.fna','w')as outfile:
				outfile.write('>'+str(mtID)+'\n'+2*mtSeq)
		except Exception as e:
			print(f'A problem occured during {self.organism_name} mtDNA acquisition!\nMore specifically:\n{str(e)}')

	def get_gDNA(self,assembly_summary:pd.DataFrame)->None:
		try:
			if os.path.exists(f'{self.out_path}mtDNA.fna')==False:
				print(f'No mtDNA sequence was found for {self.organism_name}!')
			elif os.path.exists(f'{self.out_path}gDNA.fna')==False:
				#get latest assembly
				ftp_path=assembly_summary.loc[
					assembly_summary['organism_name']==self.organism_name.replace('_',' ').capitalize()
				]['ftp_path'].tolist()[0]
				current_assembly=ftp_path.split('/')[-1]
				#download latest assembly genome
				call(
						f'wget --output-document={self.out_path}gDNA.fna.gz {ftp_path}/{current_assembly}_genomic.fna.gz -q',
						shell=True
					)
				#decompress gDNA
				call(f'gzip -d {self.out_path}gDNA.fna.gz',shell=True)
		except Exception as e:
			print(f'A problem occured during {self.organism_name} gDNA acquisition!\nMore specifically:\n{str(e)}')

	def LASTALignment(self)->None:
		try:
			if (os.path.getsize(f'{self.out_path}mtDNA.fna')>100):
				#generate LASTAL db
				call(
						f'lastdb {self.out_path}db {self.out_path}gDNA.fna',
						shell=True
					)
				#align gDNA with mtDNA
				call(
						f"""lastal -r{self.settings['match_score']} -q{self.settings['mismatch_score']} -a{self.settings['gap_opening_score']} -b{self.settings['gap_extension_score']} {self.out_path}db {self.out_path}dmtDNA.fna  > {self.out_path}aligned_dmtDNA.afa""",
						shell=True
					)
			else:
				print(f'A problem occured during {self.organism_name} db building or alignment!')
		except Exception as e:
			print(f'A problem occured during {self.organism_name} db building or alignment!\nMore specifically:\n{str(e)}')

	def process_alignment(self)->None:
		try:
			df_input=[]
			#get data from alignment file
			with open(f'{self.out_path}aligned_dmtDNA.afa')as infile:
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
			mtRecord=SeqIO.read(f"{self.out_path}mtDNA.fna", "fasta")
			mtSize=len(str(mtRecord.seq))
			size_fil=alignments['mitochondrial_start']<mtSize
			#get current assembly for filename
			#get latest assembly
			ftp_path=self.assembly_summary.loc[
				self.assembly_summary['organism_name']==self.organism_name.replace('_',' ').capitalize()
			]['ftp_path'].tolist()[0]
			current_assembly=ftp_path.split('/')[-1]
			#apply filters and write output into results folder
			alignments=alignments[size_fil]
			alignments.to_csv(f'{self.out_path}{self.organism_name}_{current_assembly}_numts.csv',header=True)
			#remove unnecessary files
			files=pd.Series(os.listdir(self.out_path))
			files[~files.str.contains("mitochondrion|assembly_summary|NUMT_finder|settings|organism|numts")].apply(lambda file: os.remove(f"{self.out_path}{file}"))
		except Exception as e:
			print(f'A problem occured during {self.organism_name} alignment processing!\nMore specifically:\n{str(e)}')