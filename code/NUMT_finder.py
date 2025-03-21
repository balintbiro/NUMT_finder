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

from NUMT_finder_class import NUMT_finder

logging.basicConfig(level=logging.INFO)
logger=logging.getLogger(__name__)

#initialize the argument parser object
parser=argparse.ArgumentParser(
		formatter_class=RawTextHelpFormatter,
		epilog="""
		Description:
		CLI application for mining NUclear MiTochondrial sequences (NUMTs) from mammalian reference genomes.

		Example usage:
		python code/NUMT_finder.py --org-path data/organism_names.txt --out-path data/
		"""
	)

parser.add_argument("--org-path",required=True,help="Path to the txt file containing organism names")
parser.add_argument("--out-path",required=True,help="Path to the output files")

#parse arguments
args=parser.parse_args()

#setup log environment
log_file=f"{args.out_path}NUMT_finder.log"
file_handler=logging.FileHandler(log_file)
file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logger.addHandler(file_handler)

logger.info("Starting NUMT mining")
start_time=time.time()

#get organism names
logger.info("Reading organism names")
organism_names=np.loadtxt(args.org_path,dtype=str,delimiter='\n')
if organism_names.size==1:
	logger.info(f"Mining is running for 1 organism")
else:
	logger.info(f"Mining is running for {len(organism_names)} organisms")

#get assembly summary file
if os.path.exists(f'{args.out_path}assembly_summary_refseq.txt')==False:
	logger.info(f"Downloading assembly summary to {args.out_path}assembly_summary_refseq.txt")
	ftp_site='https://ftp.ncbi.nlm.nih.gov/genomes/refseq/'
	filename='assembly_summary_refseq.txt'
	call(f'wget --directory-prefix={args.out_path} {ftp_site}{filename} -q',shell=True)
	logger.info(f"{args.out_path}assembly_summary_refseq.txt has been successfully downloaded")

#read assembly summary file
assembly_summary=pd.read_csv(
		f'{args.out_path}assembly_summary_refseq.txt',
		usecols=["organism_name","ftp_path","assembly_level"],
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
with open('settings.json')as infile:
	lastal_settings=infile.read()
lastal_settings=json.loads(lastal_settings)

#download mitochondrial file
if os.path.exists(f'{args.out_path}mitochondrion.1.1.genomic.fna')==False:
	logger.info(f"Downloading and uncompressing mitochondrial genome file to {args.out_path}mitochondrion.1.1.genomic.fna")
	call(
			f'wget --directory-prefix={args.out_path} https://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.1.genomic.fna.gz -q',
			shell=True
		)
	call(
			f'gzip -d {args.out_path}mitochondrion.1.1.genomic.fna.gz',
			shell=True
		)
	logger.info(f"{args.out_path}mitochondrion.1.1.genomic.fna has been successfully downloaded and uncompressed")

def main()->None:
	if organism_names.size>1:
		for organism_name in organism_names:
			NUMT_class=NUMT_finder(organism_name=organism_name,out_path=args.out_path,settings=lastal_settings,assembly_summary=assembly_summary)
			logger.info(f"mtDNA is being downloaded for {organism_name}")
			NUMT_class.get_mtDNA()
			logger.info(f"mtDNA was successfully downloaded for {organism_name}")
			logger.info(f"gDNA is being downloaded for {organism_name}")
			NUMT_class.get_gDNA(assembly_summary=assembly_summary)
			logger.info(f"gDNA was successfully downloaded for {organism_name}")
			logger.info("LAST alignment is in progress")
			NUMT_class.LASTALignment()
			logger.info(f"LAST alignment has been successfully finished for {organism_name}")
			NUMT_class.process_alignment()
	else:
		NUMT_class=NUMT_finder(organism_name=organism_names,out_path=args.out_path,settings=lastal_settings,assembly_summary=assembly_summary)
		logger.info(f"mtDNA is being downloaded for {organism_names}")
		NUMT_class.get_mtDNA()
		logger.info(f"mtDNA was successfully downloaded for {organism_names}")
		logger.info(f"gDNA is being downloaded for {organism_names}")
		NUMT_class.get_gDNA(assembly_summary=assembly_summary)
		logger.info(f"gDNA was successfully downloaded for {organism_names}")
		logger.info("LAST alignment is in progress")
		NUMT_class.LASTALignment()
		logger.info(f"LAST alignment has been successfully finished for {organism_names}")
		NUMT_class.process_alignment()
	end_time=time.time()
	logger.info(f"Pipeline completed in {end_time-start_time:.2f} seconds")
	logger.info(f"Results saved to {args.out_path}")

if __name__=="__main__":
	main()