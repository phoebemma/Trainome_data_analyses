#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 10:19:04 2022

@author: chidimma
"""

import argparse
from argparse import ArgumentParser
import os
import subprocess
from colorama import Fore, Back, Style
import fnmatch
import sys
import glob


#The following lines define the command line arguements. Two arguments 
#(input directory and output directory)will be required 
#while the rest will be optional. The user defines which analyses to run 
#and what tools to use. 
parser = argparse.ArgumentParser(prog="Trainome_Pipeline", 
description = "The Trainome_pipeline is designed to ease RNA-Seq analyses by combining different analyses tools and options .\
It is intended to make  RNA sequence analyses easier both for those with \
robust bioinformatics skills and those relatively new to bioinformatics")
Required=parser.add_argument_group("Required", "These inputs are required. That is, not optional")
Required.add_argument("-i", "--input_dir", required=True, dest="input_dir", \
help="the path to the input directory where RNA sequence data are stored in fastq.gz format")

Required.add_argument ("-o", "--output_dir", dest= "output_dir", required=True, \
help ="the path to the directory where results outputs should be saved")

parser.add_argument("-t", "--threads", dest="threads",type = int, default = 1,
help="the number of threads to be used in executing the analyses. This defaults to\
   1 if none is supplied")
parser.add_argument("-g", "--trim_galore", dest="trim_galore_trim", 
help="requests trimming using TrimGalore")

#Create a mutually exclusive option for STAR-mapping
STAR = parser.add_mutually_exclusive_group()
STAR.add_argument("-p", "--PE_STAR_mapping", dest="PE_STAR_mapping", nargs=2,
help="requests paired_end mapping using STAR: This is set to generate Bam files,\
gene as well as transcript counts. Add two paths; that to your STAR genome \
index, and the other to your primary assembly annotation in gtf format (IN THAT ORDER)")
STAR.add_argument("-q", "--SE_STAR_mapping", dest="SE_STAR_mapping", nargs=2,
help="requests single end mapping using STAR: This is set to generate Bam files,\
gene as well as transcript counts. Add two paths; that to your STAR genome \
index, and the other to your primary assembly annotation in gtf format (IN THAT ORDER)")

parser.add_argument("-a", "--alternative_splicing", dest="alternative_splicing",
nargs=2, help="requests alternative splicing analyses using Spladder. This \
analyses the presence of the following alternative splicing events in input \
data; exon skips, intron retention, alternative 3’ splice sites, alternative \
5’ splice sites,multiple (coordinated) exons skips, mutually exclusive exons ")

parser.add_argument("-k", "--kallisto", dest="kallisto_mapping",nargs=2, 
help="requests mapping using kallisto. This generates transcripts counts")

parser.add_argument("-e", "--splicing_efficiency", dest="splicing_efficiency", 
help="requests analyses of splicing efficiency using SPLICE-q. This ranks intron \
retention from 0-1. 1 being completely spliced intron and 0 complete intron retention")

parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.01", 
help="prints out the pipeline version")

#Create a mutually exculsive option for trimmomatic. Users can either perform
# paired_end trimming, or single end trimming
trimmomatic = parser.add_mutually_exclusive_group()
trimmomatic.add_argument("-r", "--PE_trimmomatic", dest="PE_trimmomatic", 
nargs=2, help = "requests trimming using Trimmomatic (Paired-end mode). Requires two \
directory paths, One to the trimmomatic jar, the other to the adapters.IN THAT ORDER\
This argument is mutually exclusive with the SE mode (-s)")
trimmomatic.add_argument("-s", "--SE_trimmomatic", dest="SE_trimmomatic", 
nargs=2, help = "requests trimming using Trimmomatic (single_end mode). Requires two \
directory paths, One to the trimmomatic jar, the other to the adapters \
(IN THAT ORDER). This argument is mutually exclusive with the PE mode (-p)")
args = parser.parse_args()

#to check if input directory exists
if not os.path.exists(args.input_dir):
 sys.exit(Fore.BLUE + "The input directory you specified does not exist or is inaccesible")

#check output directory
if not os.path.exists(args.output_dir):
 sys.exit(Fore.BLUE + "The output directory you specified does not exist or is inaccesible")

output_dir = args.output_dir
input_dir = args.input_dir
threads = args.threads
PE_STAR_mapping = args.PE_STAR_mapping
SE_STAR_mapping = args.SE_STAR_mapping
PE_trimmomatic=args.PE_trimmomatic
SE_trimmomatic=args.SE_trimmomatic
trim_galore= args.trim_galore_trim
alt_splicing = args.alternative_splicing
kallisto_quant = args.kallisto_mapping
splicing_efficiency = args.splicing_efficiency

#Grab all fastq files in the input directory
Fasta_files = []
directory = fnmatch.filter(os.listdir(input_dir), "*.fastq.gz")
for file in directory:
    Fasta_files.append(os.path.join(input_dir, file))


#Grab all forward sequences
Forward_fasta_files = []

#Grab all reverse sequences
Reverse_fasta_files = []
#for file in Fasta_files:
   # if fnmatch.fnmatch(file, "*R1_001.fastq.gz"):
       # Forward_fasta_files.append(file)
    #elif:
        #if fnmatch.fnmatch(file, "*R2_001.fastq.gz"):
            #Reverse_fasta_files.append(file)
for filename in glob.glob(os.path.join(input_dir, "*_R1_*.fastq.gz")):
    filename_match = filename.replace("_R1_", "_R2_")
    
    if not os.path.exists(filename_match):
        print("match not found for %s..." % filename)
        continue
    Forward_fasta_files.append(filename)
    Reverse_fasta_files.append(filename_match)
    
    
    
    
def fastqc(file):
    """
    This function makes a new folder called "fastqc results" (if none exists)
    in the output directory, runs quality check on the input sequences using
    fastqc and saves them in the new folder
    """
    fastq = os.path.join(output_dir, "fastqc_results")
    if not os.path.exists(fastq):
        os.makedirs(fastq)
    command= "fastqc {} -o {} -t {}".format(file,fastq, threads)
    df = subprocess.check_call(command, shell=True)
    return(df)
    print("Quality check completed, check results in %s!" %fastq)
    

def PE_trim_trimmomatic(file):
    """
    This function creates a new folder called "trimmed with trimmomatic 
    sequences" in the output directory and trims the sequences using trimmomatic 
    in paired_end mode
    """
    global trim
    trim= os.path.join(output_dir, "trimmed_with_trimmomatic_sequences")
    if not os.path.exists(trim):
        os.makedirs(trim)
    file_basename=os.path.basename(file)
    trimmomatic_jar = PE_trimmomatic[0]
    adapter_seqs = PE_trimmomatic[1]
    command ="java -jar {} PE -threads {} -basein {} -baseout {}/{} ILLUMINACLIP:{}:2:30:15:8:true MAXINFO:40:0.7 MINLEN:10 LEADING:3 TRAILING:3".format(trimmomatic_jar, threads, file,  trim, file_basename, adapter_seqs)
    df = subprocess.check_call(command, shell = True)
    return df

def SE_trim_trimmomatic(file):
    """
    This function creates a new folder called "trimmed with trimmomatic 
    sequences" in the output directory and trims the sequences using trimmomatic 
    in single_end mode
    """
    global trim_SE
    trim_SE= os.path.join(output_dir, "trimmed_with_trimmomatic_sequences")
    if not os.path.exists(trim_SE):
        os.makedirs(trim_SE)
    file_basename=os.path.basename(file)
    trimmomatic_jar = SE_trimmomatic[0]
    adapter_seqs = SE_trimmomatic[1]
    command ="java -jar {} SE -threads {} {} {}/{} ILLUMINACLIP:{}:2:30:15:8:true MAXINFO:40:0.7 MINLEN:10 LEADING:3 TRAILING:3".format(trimmomatic_jar, threads, file,  trim_SE, file_basename, adapter_seqs)
    df = subprocess.check_call(command, shell = True)
    return df


def fastqc_after_trimmomatic(file):
    """
    This function runs a quality check on trimmed sequences after trimming
    """
    fastqc_trimmed = os.path.join(output_dir, "fastqc_trimmed_sequences")
    if not os.path.exists(fastqc_trimmed):
        os.makedirs(fastqc_trimmed)
    command= "fastqc {} -o {} -t {}".format(file,fastqc_trimmed, threads)
    df = subprocess.check_call(command, shell=True)
    return(df)

    

def PE_mapping_STAR(file1, file2):
    global PE_STAR
    PE_STAR = os.path.join(output_dir, "PE_STAR_output")
    if not os.path.exists(PE_STAR):
        os.makedirs(PE_STAR)
    genome_directory= args.PE_STAR_mapping[0]
    annotation = args.PE_STAR_mapping[1]
    file_basename = os.path.basename(file1)
    command = "STAR --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts  --runDirPerm All_RWX --runThreadN {} --sjdbGTFfile {} --readFilesCommand zcat --winAnchorMultimapNmax 150 --outFilterMultimapNmax 80 --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMstrandField intronMotif --outWigType wiggle --outFilterMismatchNoverLmax 0.1 --genomeDir {} --readFilesIn {} {} --outFileNamePrefix {}/{}".format(threads,annotation, genome_directory, file1, file2, PE_STAR, file_basename)
    df = subprocess.check_call(command, shell = True)
    return(df)


def SE_mapping_STAR(file):
    global SE_STAR
    SE_STAR = os.path.join(output_dir, "SE_STAR_output")
    if not os.path.exists(SE_STAR):
        os.makedirs(SE_STAR)
    genome_directory= args.SE_STAR_mapping[0]
    annotation = args.SE_STAR_mapping[1]
    file_basename = os.path.basename(file)
    command = "STAR --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts  --runDirPerm All_RWX --runThreadN {} --sjdbGTFfile {} --readFilesCommand zcat --winAnchorMultimapNmax 150 --outFilterMultimapNmax 80 --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMstrandField intronMotif --outWigType wiggle --outFilterMismatchNoverLmax 0.1 --genomeDir {} --readFilesIn {} --outFileNamePrefix {}/{}".format(threads,annotation, genome_directory, file, SE_STAR, file_basename)
    df = subprocess.check_call(command, shell = True)
    return(df)



print("........................................................................")
print("Checking the quality of your samples using fastqc")

for file in Fasta_files:
    fastqc(file)
    continue

if PE_trimmomatic:
    for file in Fasta_files:
        PE_trim_trimmomatic(file)

elif SE_trimmomatic:
    for file in Fasta_files:
        SE_trim_trimmomatic(file)
        
        continue
    
#creating a list of the trimmed sequences using trimmomatic in PE mode
trimmed_files_PE= []
trimmed_1P = []
trimmed_2P =[]
for trimmed in trim:
    trimmed_files_PE.append(trimmed)
for trimmed in glob.glob(os.path.join(trim, "*_1P.fastq.gz")):
   trimmed_match = trimmed.replace("_1P", "_2P")
   if not os.path.exists(trimmed_match):
        print("match not found for %s..." % trimmed)
        continue
   trimmed_1P.append(trimmed)
   trimmed_2P.append(trimmed_match)

#creating a list of the trimmed sequences using trimmomatic in the SE mode
trimmed_files_SE=[]
for trimmed in trim_SE:
    trimmed_files_SE.append(trimmed)
#rUNNING A QUALITY CHECK OF THE TRIMOMATIC OUTPUTS  
if PE_trimmomatic:
    for file in trimmed_files_PE:
        fastqc_after_trimmomatic()
  
elif SE_trimmomatic:
    for file in trimmed_files_SE:
        fastqc_after_trimmomatic()
        continue

# checks if the user requested for trimming and STAR mapping
#if user requests for both, the input sequences for the mapping
# are the trimmed sequences. If user requests only for mapping , input sequences
#will be the raw fastq files
if PE_trimmomatic and PE_STAR_mapping:
    #use the trimmed sequences as input directory
    for file, file2 in zip(trimmed_1P, trimmed_2P):
        PE_mapping_STAR(file, file2)
else:
    if PE_STAR_mapping and not PE_trimmomatic:
        for file, file2 in zip(Forward_fasta_files, Reverse_fasta_files):
            PE_mapping_STAR(file, file2)
            continue

# checks for single end trimming and mapping as above
if SE_trimmomatic and SE_STAR_mapping:
        for file in trimmed_files_SE:
            SE_mapping_STAR(file)
else:
    if SE_STAR_mapping and not SE_trimmomatic:
        for file in Fasta_files:
            SE_mapping_STAR()



