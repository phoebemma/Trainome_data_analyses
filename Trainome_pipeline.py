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

parser = argparse.ArgumentParser(prog="Trainome_Pipeline", description = "The Trainome_pipeline is designed to ease RNA-Seq analyses by combining different analyses tools and options . It is intended to make  RNA sequence analyses easier both for those with robust bioinformatics skills and those relatively new to bioinformatics")
parser.add_argument("-i", "--input_dir", dest="input_dir", required=True, help="the path to the input directory where RNA sequence data are stored in fastq.gz")
parser.add_argument("-o", "--output_dir", dest="output_dir", required=True, help="the path to the directory where results outputs should be saved")
parser.add_argument("-t", "--threads", dest="threads",type=int, required=True, help="the number of threads to be used in executing the analyses")
parser.add_argument("-r", "--trimmomatic",action="store_true", dest="trimmomatic_trim", help="requests trimming using Trimmomatic")
parser.add_argument("-g", "--trim_galore", dest="trim_galore_trim", help="requests trimming using TrimGalore")
parser.add_argument("-s", "--STAR_mapping", dest="STAR_mapping", help="requests mapping using STAR: This is set to generate Bam files, gene as well as transcript counts")
parser.add_argument("-a", "--alternative_splicing", dest="alternative_splicing", help="requests alternative splicing analyses using Spladder. This analyses the presence of the following alternative splicing events in input data; exon skips, intron retention, alternative 3’ splice sites, alternative 5’ splice sites,multiple (coordinated) exons skips, mutually exclusive exons ")
parser.add_argument("-k", "--kallisto", dest="kallisto_mapping", help="requests mapping using kallisto. This generates transcripts counts")
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
STAR_mapping = args.STAR_mapping
trimmomatic = args.trimmomatic_trim
trim_path = "/datastore/Chidimma/tools/Trimmomatic-0.39"
fastqc = "fastqc $file -o$output -t $threads"
if args.trimmomatic_trim:
    fastq= os.path.join(output_dir, "fastqc_results")
  ##     os.mkdir(fastq)
    #for file in os.listdir(input_dir):
     #   if fnmatch.fnmatch(file, "*.fastq.gz"):
      #      subprocess.call(["sh", "./trim.sh"])
#else:
 #   pass
#if args.STAR_mapping
#os.environ['threads']=threads
#os.environ['output']=output
#Grab all fastq files in the input directory
Fasta_files = []
#for file in fnmatch.filter(os.listdir(input_dir), "*.fastq.gz"):
directory = fnmatch.filter(os.listdir(input_dir), "*.fastq.gz")
for file in directory:
    Fasta_files.append(os.path.join(input_dir, file))
#env = {"file": "file", "threads" : "threads", "fastq" : "fastq"}
def fastqc(file):
    command= "fastqc {} -o {} -t {}".format(file, output_dir, threads)
    df = subprocess.check_call(command, shell=True)
    return(df)

if trimmomatic:
    fastq = os.path.join(output_dir, "fastqc_results")
    if not fastq:
        os.mkdir(fastq)    

    for file in  Fasta_files:
        fastqc(file)