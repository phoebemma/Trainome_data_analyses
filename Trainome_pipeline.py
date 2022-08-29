#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 10:19:04 2022

@author: chidimma
"""

import argparse
from argparse import ArgumentParser
import os
from os import path
import subprocess
from colorama import Fore, Back, Style
import pathlib
import fnmatch
import sys
import re

parser = argparse.ArgumentParser(prog="Trainome_Pipeline", description = "The Trainome_pipeline is designed to ease RNA-Seq analyses by combining different analyses tools and options . It is intended to make  RNA sequence analyses easier both for those with robust bioinformatics skills and those relatively new to bioinformatics")
parser.add_argument("-i", "--input_dir", dest="input_dir", required=True, help="the path to the input directory where RNA sequence data are stored in fastq.gz")
parser.add_argument("-o", "--output_dir", dest="output_dir", required=True, help="the path to the directory where results outputs should be saved")
parser.add_argument("-t", "--threads", dest="threads",type=int, required=True, help="the number of threads to be used in executing the analyses")
parser.add_argument("-r", "--trimmomatic", dest="trimmomatic_trim", nargs=2, help="requests trimming using Trimmomatic (Paired-end mode). Requires two directory paths, One to the trimmomatic jar, the other to the Adapters path")
parser.add_argument("-g", "--trim_galore", dest="trim_galore_trim", help="requests trimming using TrimGalore")
STAR_group = parser.add_mutually_exclusive_group()
STAR_group.add_argument("-p", "--PE_STAR_mapping", dest="PE_STAR_mapping", nargs=2, help="requests paired_end mapping using STAR: This is set to generate Bam files, gene as well as transcript counts. Add two paths; that to your STAR genome index, and the other to your primary assembly annotation in gtf format")
STAR_group.add_argument("-s", "--SE_STAR_mapping", dest="SE_STAR_mapping", nargs=2, help="requests single_end mapping using STAR: This is set to generate Bam files, gene as well as transcript counts. Add two paths; that to your STAR genome index, and the other to your primary assembly annotation in gtf format")
parser.add_argument("-a", "--alternative_splicing", dest="alternative_splicing",nargs=2, help="requests alternative splicing analyses using Spladder. This analyses the presence of the following alternative splicing events in input data; exon skips, intron retention, alternative 3’ splice sites, alternative 5’ splice sites,multiple (coordinated) exons skips, mutually exclusive exons ")
parser.add_argument("-k", "--kallisto", dest="kallisto_mapping",nargs=2, help="requests mapping using kallisto. This generates transcripts counts")
parser.add_argument("-e", "--splicing_efficiency", dest="splicing_efficiency", help="requests analyses of splicing efficiency using SPLICE-q. This ranks intron retention from 0-1. 1 being completely spliced intron and 0 complete intron retention")
parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.01", help="prints out the pipeline version")
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
trimmomatic=args.trimmomatic_trim
trim_galore= args.trim_galore_trim
alt_splicing = args.alternative_splicing
kallisto_quant = args.kallisto_mapping
splicing_efficiency = args.splicing_efficiency

#fastq=subprocess.check_call(["mkdir {}/fastqc_results".format(output_dir)], shell =True)
#fastq = os.path.join(output_dir, "fastqc_results")
#if not fastq:
    #os.makedirs(fastq)
#if args.trimmomatic_trim:
    #fastq= os.path.join(output_dir, "fastqc_results")
  ##     os.mkdir(fastq)
    #for file in os.listdir(input_dir):
     #   if fnmatch.fnmatch(file, "*.fastq.gz"):
      #      subprocess.call(["sh", "./trim.sh"])
#else:
#fastq = os.path.join(output_dir, "fastqc_results")
#if not os.path.exists(fastq):
#    os.makedirs(fastq)

#Grab all fastq files in the input directory
Fasta_files = []
#for file in fnmatch.filter(os.listdir(input_dir), "*.fastq.gz"):
directory = fnmatch.filter(os.listdir(input_dir), "*.fastq.gz")
for file in directory:
    Fasta_files.append(os.path.join(input_dir, file))


#Grab all forward sequences
Forward_fasta_files = []
Reverse_fasta_files = []
for file in Fasta_files:
    if fnmatch.fnmatch(file, "*R1_001.fastq.gz"):
        Forward_fasta_files.append(file)
    else:
        Reverse_fasta_files.append(file)
    

def fastqc(file):
    #make a folder for the results if none exists
    fastq = os.path.join(output_dir, "fastqc_results")
    if not os.path.exists(fastq):
        os.makedirs(fastq)
    #fastq= os.makedirs("fastq_results")
    command= "fastqc {} -o {} -t {}".format(file,fastq, threads)
    #command= "fastqc {} -o {} -t {}".format(file,fastq, threads)
    df = subprocess.check_call(command, shell=True)
    return(df)

def trim_trimmomatic(file):
    trim= os.path.join(output_dir, "trimmmed_with_trimmomatic_sequences")
    if not os.path.exists(trim):
        os.makedirs(trim)
    file_basename=os.path.basename(file)
    trimmomatic_jar = args.trimmomatic_trim[0]
    adapter_seqs = args.trimmomatic_trim[1]
    command ="java -jar {} PE -threads {} -basein {} -baseout {}/{} ILLUMINACLIP:{}:2:30:15:8:true MAXINFO:40:0.7 MINLEN:10 LEADING:3 TRAILING:3".format(trimmomatic_jar, threads, file,  trim, file_basename, adapter_seqs)
    df = subprocess.check_call(command, shell = True)
    return df

def fastqc_after_trimmomatic(file):
    #make fastqc_trimmed available outside the function
    global fastqc_trimmed
    fastqc_trimmed = os.path.join(output_dir, "fastqc_trimmed_sequences")
    if not os.path.exists(fastqc_trimmed):
        os.makedirs(fastqc_trimmed)
    command= "fastqc {} -o {} -t {}".format(file,fastqc_trimmed, threads)
    df = subprocess.check_call(command, shell=True)
    return(df)
    

def PE_STAR_mapping(file):
    global PE_STAR
    PE_STAR = os.path.join(output_dir, "PE_STAR_output")
    if not os.path.exists(PE_STAR):
        os.makedirs(PE_STAR)
    genome_directory= args.PE_STAR_mapping[0]
    annotation = args.PE_STAR_mapping[1]
    file_basename = os.path.basename(file)
    file1=fnmatch.fnmatch(file, "[!]")
    command = "STAR --outSAMattributes All --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts  --runDirPerm All_RWX --runThreadN {} --sjdbGTFfile {} --readFilesCommand zcat --winAnchorMultimapNmax 150 --outFilterMultimapNmax 80 --outReadsUnmapped Fastx --outMultimapperOrder Random --outWigType wiggle --outFilterMismatchNoverLmax 0.1 --genomeDir {} --readFilesIn {} {} --outFileNamePrefix {}/{}".format(threads,annotation, genome_directory, file1, file2, PE_STAR, file_basename)
    
if trimmomatic:
    for file in Fasta_files:
        #Choose only the forward sequences, this is because the trimmomatic command is using the "-basein" option that automatically determines the reverse file
        if fnmatch.fnmatch(file, "*R1_001.fastq.gz"):
            trim_trimmomatic(file)
    for file in fastqc_trimmed:
        if fnmatch.fnmatch(file, "fastq.gz"):
            fastqc_after_trimmomatic(file)
     
    #UNCOMMENT FOR FASTQC LATER
    for file in  Fasta_files:
        fastqc(file)