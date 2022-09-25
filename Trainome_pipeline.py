#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 08:21:48 2022

@author: chidimma
"""
#import configparser
import argparse, configparser
from argparse import ArgumentParser
import os
import subprocess
from colorama import Fore, Back, Style
import fnmatch
import sys
import glob

DEFAULT="trimmomatic.ini"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Trainome_Pipeline", 
    description = "The Trainome_pipeline is designed to ease RNA-Seq analyses by combining different analyses tools and options .\
    It is intended to make  RNA sequence analyses easier both for those with \
    robust bioinformatics skills and those relatively new to bioinformatics")
   

# Defining the required arguments
Required=parser.add_argument_group("Required", "These inputs are required. That is, not optional")
Required.add_argument("-i", "--input_dir", required=True,  metavar = ".",dest="input_dir",\
help="the path to the input directory where RNA sequence data are stored in fastq.gz format")

Required.add_argument ("-o", "--output_dir", dest= "output_dir",  metavar = ".", required=True,
help ="the path to the directory where results outputs should be saved") 
    
#Number of threads. Defaults to 1
parser.add_argument("-t", "--threads", dest="threads",type = int, metavar = "",  default = 1,
help="the number of threads to be used in executing the analyses. This defaults to 1 if none is supplied")

#Define the mode of experiment
parser.add_argument("mode", choices= ["PE", "SE"], 
 help = "requests for paired-end (PE) or single_end (SE) sequence analyses", )

#add  mutually exclusive arguments of for trimming with trimmomatic
trimmomatic_trim = parser.add_mutually_exclusive_group()
trimmomatic_trim.add_argument("-r", "--trimmomatic_trim", nargs =2, dest = "trimmomatic", 
help = "requests trimming using Trimmomatic. Requires two paths. The first to the  \
 trimmomatic jar, the others to the adapters.", metavar = "")

trimmomatic_trim.add_argument("--trim_commands", dest = "trim_commands", nargs = "?", type = argparse.FileType("r"))

#add  mutually exclusive arguments for STAR_mapping
Star_mapping = parser.add_mutually_exclusive_group()
Star_mapping.add_argument("-s", "--STAR_mapping", dest = "STAR_mapping", metavar = "", nargs=2, 
help = "requests sequence mapping using STAR. The defaults setting will generate\
bam files and counts of reads per gene. Two paths are needed, first path to user's STAR \
genome index, the second to the primary assembly annotation in gtf format" )

Star_mapping.add_argument("--STAR_commands",dest = "STAR_commands", nargs = "?", type = argparse.FileType("r"))

parser.add_argument("-k", "--kallisto", dest="kallisto_quant", metavar = "", nargs=1, 
help="requests mapping using kallisto. This generates transcripts counts. \
The path to a kallisto index is needed")
parser.add_argument("-e", "--splicing_efficiency", dest = "splicing_efficiency",
metavar = "", nargs = 1, help = "requests using Splice_q to analyse the splicing \
efficiency of the data. This requires one input; the path to the primary \
assembly annotation file (in gtf format)")
parser.add_argument("--Rsem", dest = "Rsem_calc_exp", metavar = "", nargs = 1,
help = "requests the use of RSEM to quantify transcripts. This requires \
    the path to the RSEM index")    

parser.add_argument("--salmon", dest= "salmon_quant", metavar = "", nargs = 1,
help = "requests the use of Salmon to quantify transcripts. Requires the path \
    to a salmon index")
    
parser.add_argument("--stringtie", dest = "stringtie", metavar = "", nargs = 1,
 help = "requests the use of Stringtie to quantify transcripts. It takes as input\
     sorted bam files . It requires the path to a genome annotation file in gtf format")
    
args = parser.parse_args()
if args.trim_commands:
    #config = configparser.ConfigParser()
    #config.read(DEFAULT)
    #defaults = {}
    #defaults.update(dict(config.items("USER_PREFERENCES")))
    #parser.set_defaults(**defaults)
    #args = parser.parse_args()
    with open(args.trim_commands, "r") as f:
        config = configparser.ConfigParser()
        config.read([args.trim_file])
    
    

#Parameters to override the args arguments
output_dir = args.output_dir
input_dir = args.input_dir
threads = args.threads
splicing_efficiency = args.splicing_efficiency
Rsem_calc_exp = args.Rsem_calc_exp
stringtie = args.stringtie

if args.mode == "PE":
    PE_trimmomatic=args.trimmomatic
    trim_commands = args.trim_commands
    PE_STAR_mapping = args.STAR_mapping
    STAR_commands = args.STAR_commands
    PE_kallisto_quant = args.kallisto_quant
    PE_salmon_quant = args.salmon_quant
   
else:
    if args.mode == "SE":
        SE_trimmomatic = args.trimmomatic
        trim_commands = args.trim_commands
        SE_STAR_mapping = args.STAR_mapping
        STAR_commands = args.STAR_commands
        SE_kallisto_quant = args.kallisto_quant
        SE_salmon_quant = args.salmon_quant
        

#Check if input and output directories exist
if not os.path.exists(input_dir):
 sys.exit(Fore.BLUE + "The input directory you specified does not exist or is inaccesible")

#check output directory
if not os.path.exists(output_dir):
 sys.exit(Fore.BLUE + "The output directory you specified does not exist or is inaccesible")


#Grab all fastq files in the input directory
Fasta_files = []
directory = fnmatch.filter(os.listdir(input_dir), "*.fastq.gz")
for file in directory:
    Fasta_files.append(os.path.join(input_dir, file))
    
    
#Create list for all forward sequences
Forward_fasta_files = []

#Create list for all reverse sequences
Reverse_fasta_files = []
#for file in Fasta_files:
   # if fnmatch.fnmatch(file, "*R1_001.fastq.gz"):
  #      Forward_fasta_files.append(file)
 #   else:
  #      if fnmatch.fnmatch(file, "*R2_001.fastq.gz"):
      #      Reverse_fasta_files.append(file)
for filename in glob.glob(os.path.join(input_dir, "*_R1_*.fastq.gz")):
    filename_match = filename.replace("_R1_", "_R2_")
    
    if not os.path.exists(filename_match):
        print("match not found for %s..." % filename)
        continue
    Forward_fasta_files.append(filename)
    Reverse_fasta_files.append(filename_match)
print(Forward_fasta_files)
print(Reverse_fasta_files)
#Define a function that runs quality check on all input files 

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
    
def samtools_indexing(file):
        if fnmatch.fnmatch(file, "*sortedByCoord.out.bam"):
            indexing = "samtools index {}".format (file)
            df = subprocess.check_call(indexing, shell = True)
            return(df)
            
    
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

def Splice_q(file):
    global Splice_file
    Splice_file = os.path.join(output_dir, "SPLICE_q_outputs")
    if not os.path.exists(Splice_file):
        os.makedirs(Splice_file)
    annotation = splicing_efficiency[0]
    file_basename = os.path.basename(file)
    indiv_folder = os.path.join(Splice_file, file_basename + ".tsv")
    command = "SPLICE-q.py -b {} -g {} -p {} -o {}". format(file, 
    annotation, threads, indiv_folder)
    df = subprocess.check_call(command, shell = True)
    return(df)

def Rsem_quant(file):
    """
    This creates a file called RSEM_outputs and runs a commandline argument that \
    quantifies the transcripts in a ".gzAligned.toTranscriptome.out.bam" file
    """
    global RSEM
    RSEM = os.path.join(output_dir, "RSEM_outputs")
    if not os.path.exists(RSEM):
        os.makedirs(RSEM)
    file_basename=os.path.basename(file)
    indiv_folder = os.path.join(RSEM, file_basename)
    annotation = Rsem_calc_exp[0]
    command = "rsem-calculate-expression --paired-end --alignments --append-names \
    -p {} {} {} {}".format(threads, file, annotation, indiv_folder)
    df = subprocess.check_call(command, shell = True)
    return(df)
    

def Stringtie_quant(file):
    global Stringtie_output
    Stringtie_output = os.path.join(output_dir, "Stringtie_outputs")
    if not os.path.exists(Stringtie_output):
        os.makedirs(Stringtie_output)
    path = "/datastore/Chidimma/tools/stringtie-2.2.1.Linux_x86_64/"
    annotation =stringtie[0]
    file_basename = os.path.basename(file)
    indiv_folder_gtf = os.path.join(Stringtie_output, file_basename + ".gtf")
    indiv_folder_tab = os.path.join(Stringtie_output, file_basename + ".tab")
    command = "{}stringtie -o {} -A {} -G {} -B {}".format(path, indiv_folder_gtf,
    indiv_folder_tab, annotation, file)
    df = subprocess.check_call(command, shell = True)
    return(df)


#Creating a class object fro paired-end sequence analyses
class Paired_end:
    def trimmomatic(file):
        """
        This function creates a new folder called "trimmed with trimmomatic 
        sequences" in the output directory and trims the sequences using trimmomatic 
        in paired_end mode
        """
        global trim
        trim= os.path.join(output_dir, "PE_trimmed_sequences")
        if not os.path.exists(trim):
            os.makedirs(trim)
        file_basename=os.path.basename(file)
        trimmomatic_jar = PE_trimmomatic[0]
        adapter_seqs = PE_trimmomatic[1]
        command ="java -jar {} PE -threads {} -basein {} -baseout\
        {}/{} ILLUMINACLIP:{}:2:30:15:8:true MAXINFO:40:0.7 MINLEN:10 \
            LEADING:3 TRAILING:3".format(trimmomatic_jar, threads, file,  
            trim, file_basename, adapter_seqs)
        df = subprocess.check_call(command, shell = True)
        return df
    
    
    
    def STAR_mapping(file1, file2):
        
        """
        This function creates a new folder in the user's output directory
        named PE_STAR_output, then uses STAR to mapp and analyse the RNA-Seq data
        . The function also uses samtools to index the "sortedByCoord.out.bam" files"""
        global PE_STAR
        PE_STAR = os.path.join(output_dir, "STAR_outputs_PE")
        if not os.path.exists(PE_STAR):
            os.makedirs(PE_STAR)
        genome_directory= PE_STAR_mapping[0]
        annotation = PE_STAR_mapping[1]
        file_basename = os.path.basename(file1)
        command = "STAR --outSAMattributes All --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM GeneCounts  \
        --runDirPerm All_RWX --runThreadN {} \
        --sjdbGTFfile {} --readFilesCommand zcat --winAnchorMultimapNmax 150 \
        --outFilterMultimapNmax 80 --outReadsUnmapped Fastx \
        --outMultimapperOrder Random --outSAMstrandField intronMotif \
        --outWigType wiggle --outFilterMismatchNoverLmax 0.1 --genomeDir {} \
        --readFilesIn {} {} --twopassMode Basic --outFileNamePrefix {}/{}".format(threads, annotation,\
        genome_directory, file1, file2, PE_STAR, file_basename)
        df = subprocess.check_call(command, shell = True)
        return(df)
    
    
    def kallisto(file1, file2):
        global kallisto_outputs
        kallisto_outputs = os.path.join(output_dir, "kallisto_outputs")
        if not os.path.exists(kallisto_outputs):
            os.makedirs(kallisto_outputs)
        file_basename = os.path.basename(file1)
        indiv_folder = os.path.join(kallisto_outputs, file_basename + "_dir")
        if not os.path.exists(indiv_folder):
            os.makedirs(indiv_folder)
        index =  PE_kallisto_quant[0]
        command = "kallisto quant -i {} -o {} {} {} -t {} --plaintext".format(index, 
        indiv_folder, file1, file2, threads)
        df= subprocess.check_call(command, shell =True)
        os.chdir(indiv_folder)
        os.renames("abundance.tsv", file_basename + ".tsv")
       
        return(df)
        
    def salmon(file1, file2):
        global PE_Salmon_outputs
        PE_Salmon_outputs = os.path.join(output_dir, "Salmon_outputs_PE")
        if not os.path.exists(PE_Salmon_outputs):
            os.makedirs(PE_Salmon_outputs)
        path = "/datastore/Chidimma/salmon-1.9.0_linux_x86_64/bin/" 
        index = PE_salmon_quant[0]
        file_basename=os.path.basename(file1)
        indiv_folder = os.path.join(PE_Salmon_outputs, file_basename + "_dir")
        if not os.path.exists(indiv_folder):
            os.makedirs(indiv_folder)
        command = "{}salmon quant -i {} -l A --writeUnmappedNames -1 {} -2 {} -z \
        -p {} -o {} ". format(path, index, file1, file2, threads, indiv_folder)
        df = subprocess.check_call(command, shell = True)
        os.chdir(indiv_folder)
        os.renames("quant.sf", file_basename + ".sf")
   
        return(df)


#First step once a folder containing fasta files is given as input
for file in Fasta_files:
    fastqc(file)
    continue

if PE_trimmomatic:
    for file in Fasta_files:
        #Choose only the forward sequences, this is because the trimmomatic command is using the "-basein" option that automatically determines the reverse file
        if fnmatch.fnmatch(file, "*R1_001.fastq.gz"):
            Paired_end.trimmomatic(file)
            
            #creating a list of the trimmed sequences using trimmomatic in PE mode

            trimmed_files_PE= []
            trimmed_1P = []
            trimmed_2P =[]
            #for trimmed in trim:
                #if fnmatch.fnmatch(trimmed, "*.fastq.gz"):
                    #trimmed_files_PE.append(trimmed)
            for trimmed in glob.glob(os.path.join(trim, "*.fastq.gz")):
                trimmed_files_PE.append(trimmed)
           
            # creating lists for both forward and reverse sequences
            for trimmed in glob.glob(os.path.join(trim, "*_1P.fastq.gz")):
               trimmed_match = trimmed.replace("_1P", "_2P")
               if not os.path.exists(trimmed_match):
                    print("match not found for %s..." % trimmed)
                    continue
               trimmed_1P.append(trimmed)
               trimmed_2P.append(trimmed_match)
               print(trimmed_1P)
               print(trimmed_2P)
               print(trimmed_files_PE)
               continue
        

    #Quality check on the trimmed sequences
if PE_trimmomatic:
    for file in trimmed_files_PE:
        if fnmatch.fnmatch(file, "*.fastq.gz"):
            fastqc_after_trimmomatic(file)
            continue
    
    
if PE_trimmomatic and PE_STAR_mapping:
    if not os.path.exists(PE_STAR_mapping[0]):
        sys.exit(Fore.BLUE + "Your STAR index does not exist. Check the given path")

    #To check if the primary assemble annotation file exists and is in gtf  format
    gtf=os.path.exists(PE_STAR_mapping[1])
    if not os.path.exists(PE_STAR_mapping[1]) and not gtf.endswith(".gtf"):
        sys.exit(Fore.BLUE + "Your annotation file path does not exist \
        or is not in the correct format (gtf)")
    #use the trimmed sequences as input directory
    for x, y in zip(trimmed_1P, trimmed_2P):
        Paired_end.STAR_mapping(x, y)
else:        
    if PE_STAR_mapping and not PE_trimmomatic:
        for x, y in zip(Forward_fasta_files, Reverse_fasta_files):
            Paired_end.STAR_mapping(x, y)
            continue
            
            #Creating a list for the STAR outputs
STAR_sortedByCoord = []
STAR_AlignedToTranscriptome= []
for bamfiles in glob.glob(os.path.join(PE_STAR, "*sortedByCoord.out.bam")):
    STAR_sortedByCoord.append(bamfiles)
           
for bamfile in glob.glob(os.path.join(PE_STAR, "*Aligned.toTranscriptome.out.bam")):
    STAR_AlignedToTranscriptome.append(bamfile)
    continue
print(STAR_sortedByCoord)
print(STAR_AlignedToTranscriptome)
        
#Indexing the "sortedByCoord .bam" files     
for file in STAR_sortedByCoord:
   samtools_indexing(file)
   continue



if PE_kallisto_quant and PE_trimmomatic:
    if not os.path.exists(PE_kallisto_quant[0]):
        sys.exit(Fore.RED + "Your kallisto index does not exist. Check the given path")
               
    for x, y in zip(trimmed_1P, trimmed_2P):
        Paired_end.kallisto(x, y)
        
        #Make a list containing the tsv files from kallisto
      
        
        
else:
    if PE_kallisto_quant and not PE_trimmomatic:
        for x, y in zip(Forward_fasta_files, Reverse_fasta_files):
            Paired_end.kallisto(x, y)
            continue

#kallisto_PE_outputs = []
#or dirs, files in os.walk(kallisto_outputs):
   #   for file in files:
  #        if file.endswith(".tsv"):
  #            kallisto_PE_outputs.append(file)
#print(kallisto_PE_outputs)
     

            
        
#Run SPLICE_q if requested
if splicing_efficiency and PE_STAR_mapping:
    for file in STAR_sortedByCoord:
        Splice_q(file)
        
else:
    if splicing_efficiency and  not PE_STAR_mapping:
        sys.exit(Fore.BLUE + "Splicing efficiency needs bam files as input. Please perform STAR mapping")


#quantify using RSEM if demanded.
if Rsem_calc_exp:
    for file in STAR_AlignedToTranscriptome:
         Rsem_quant(file)
         continue
     
if PE_salmon_quant:
    if not os.path.exists(PE_salmon_quant[0]):
        sys.exit(Fore.RED + "Your salmon index does not exist. Check the given path")
               
    for x, y in zip(trimmed_1P, trimmed_2P):
        Paired_end.salmon(x, y)
        continue


#quantify using Stringtie
if stringtie:
    for file in STAR_sortedByCoord:
        Stringtie_quant(file)
        continue
        
    