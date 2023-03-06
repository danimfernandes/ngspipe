#!/usr/bin/python3

__author__ = 'daniel'

## Changelog ##
# v0.871
# - Added hapConX() to estimate contamination with hapConX
#
# v0.870
# - Updated sexing() to Python3
#
# v0.869
# - Fixed bug in bamDemux()
#
# v0.868
# - Fixed bamDemux() index order and running details
# - Fixed end_analysis() terminal deamination reporting
#
# v0.867
# - Added bamDemux() to demultiplex Illumina reads stored as a BAM file
#
# v0.865
# - Added schmutzi() to estimate contamination with schmutzi
# - Fixed close() for hirisplex embedded list
#
# v0.864
# - Added pileupCaller() to call pseudo-haploid genotypes, and replace pca_lc()
# - Added pileup2ped() to convert samtools mpileup data from pileupCaller() into PLINK format
#
# v0.863
# - Added changeBamChrNotation() to change BAM chromosome notation from 1,2,3 to chr1,chr2,chr3, and vice-versa
#
# v0.862
# - Added choice to use Picardtools in batchBamToFastq
# - Fixed coverage() when MT coverage was 0
#
# v0.861
# - Fixed some bugs generatng ped and map files
#
# v0.860
# - Fixed some calico() bugs, allowed multipleLanes=False
#
# v0.859
# - Added batchBAMtoFASTQ() to batch convert BAM to FASTQ
# - Added imputeWithBeagle() for imputation and phasing
#
# v0.858
# - Added yleaf() for Y haplogroup assessment
#
# v0.857
# - Updated GATK callers to version 4.1.2.0
# - Updated PICARD callers to version 2.20.2
#
# v0.856
# - Added SAMTOOLS and GATK options to makeindexgenome()
# - Corrected all occurrences that used quality thresholds to accept changes from v.0.855
# - Grouped functions by function at the end of the script
#
# v0.855
# - Added multi-threading in bwaAlign() as an argument
# - Added quality score filter for tagAndConvertSb() as an argument
# - Incorporated option to choose samtools version (1.X.X or 0.1.XX) in samToolsSort()
# - Removed plenty of old #commented code
#
# v0.854
# - Faster fastqdemux() with mismatches and dual indexing allowed
#
# v0.853
# - Less buggy Yfitter() code
#
# v0.852
# - Added fastqdemux() to demultiplex Illumina sequencing data in *.fastq[.gz] format
# - Added mergedEnd_Analysis() to create a Summary file from merged BAM files
#
# v0.851
# - Added option to use cutadapt() directly from compressed fastqs (*.gz)
# - Some calico() integration bugs corrected
#
# v0.850
# - Updated chrConcat()
#
# v0.841
# - Added calico() to estimate mitochondrial contamination
#
# v0.832
# - Added mergeBamsBatch() to batch merge BAM files from multi-lane sequencing runs
#
# v0.820
# - Added new function samToolsSort131() to be used with samtools v.1.3.1
#
# v0.790
# - Updated code to create *.map files in pca_lc(), up to 50 times faster
#
# v0.780
# - Added angsd_contamination() to integrate ANGSD to estimate X-chr contamination in males
#
# v0.770
# - Added editHeaderBAMbatch() to edit headers and read group information for all BAM files in folder
#
# v0.760
# - Added reorderBAM() to reorder contigs in BAM files that do not match the order on the desired reference
#
# v0.750
# - Added compatibility with paired end sequencing data for cutadapt() and tagAndConvertSb()
#
# v0.700
# - Added trimFastq() to trim bases at the ends of reads using SEQTK]
# - Added rescaleMapDamage() to rescale base quality scores with MAPDAMAGE algorithm


import os
import sys
import shutil
import time
import random
import re
import inspect
import regex
import gzip


# start_time = time.time()


# indir = '/home/dfernandes/Google Drive/Scripts/NGS'
# print('PATH.DIRNAME - ',os.path.dirname(indir)) ## Dirname gives mother directory of 'path'. If 'path' is a folder, gives mother folder of that folder
#
#
# cwd = os.getcwd()
# print('CWD - ',cwd)
#
# print('LISTDIR - ',cwd,os.listdir(cwd))
#
# newdir = cwd+'/Test'
# os.mkdir(newdir)
# print('MKDIR = CWD?',newdir+'?',newdir == cwd)
#
# os.chdir(newdir)
# cwd = os.getcwd()
# print('CHDIR = CWD?',newdir+'?',os.getcwd() == cwd)
# print('LISTDIR - ',cwd,os.listdir(os.getcwd()))
#
# print('DIRNAME - ',os.path.dirname((cwd)))
# os.chdir(os.path.dirname(cwd))
# print('CWD -',os.getcwd())
# cwd = os.getcwd()
# print('LISTDIR - ',cwd,os.listdir(cwd))
# print('WALK - ',cwd,os.walk(cwd))
#
# os.renames('Test','TestRe')
# print('RENAME - ',os.listdir(cwd))
# 1
# input('Delete '+newdir+'?')
# os.rmdir(cwd+'/TestRe')
#
# print(os.getcwd())

def output(filename, seqPlat, seqTypeEnd):
    """
    Function to be used as parameter section for the analysis. Must always run.

    :param filename:
    (str)   Name of file to dump run info and commands

    :param seqPlat:
    (str)   Sequencing platform. Accepts "MiSeq" and "NextSeq"

    :param seqTypeEnd:
    (str)   Sequencing type. Accepts "singleEnd" and "pairedEnd"
    :return:
    """
    global fout
    global seqPlatform
    global seqType
    seqType = seqTypeEnd
    seqPlatform = seqPlat
    cwd = os.getcwd()
    global originalWD
    originalWD = cwd
    outputcount = 0
    global mapDamageRescaleSim
    global seqtkSim
    mapDamageRescaleSim = 0
    seqtkSim = 0
    for i in os.listdir(cwd):
        if re.search(".mapdamaged.bam$", i) != None:
            mapDamageRescaleSim = 1
        if re.search(".seqtk.fastq$", i) != None:
            seqtkSim = 1
        if "_commands" in i and ".txt" in i:
            outputcount += 1
    if outputcount >= 1:
        fout = open(filename + '_commands' + str(outputcount) + '.txt', 'w')

    else:
        fout = open(filename + '_commands.txt', 'w')

    print("seqtkSim == " + str(seqtkSim))
    print("mapDamageRescaleSim == " + str(mapDamageRescaleSim))

    fout.write('>>> ' + time.strftime(
        "%d/%m/%Y\t%H:%M:%S %Z") + "\n>>> User: Daniel Fernandes" + "\n>>> Description: Listing of commands for NGS-" + seqPlatform + " analysis\n")
    # for item in samples_list:
    #     fout.write('\t'+item+'\n')
    fout.write("\n")


def fastqcAnalysis():
    # def fastqcAnalysis(moveGzFiles = False):
    """
    TODO - Port to FastQC versions > 0.10.1

    ([boolean]) -> folder

    For current folder, executes FastQC analysis for each 'fastq.gz' file, creating individual folders.
    Second argument 'moveGzFiles' is optional and is 'False' by default.
    Tested and working with FastQC version 0.10.1. Newer versions such as 0.11.3 do not work yet in providing the info for the "End_analysis".

    Input:
    foo_S11_L001_R1_001.fastq.gz

    Output:
    foo_S11_L001_R1_001.fastqc.html

    Example:
    >>> fastqcAnalysis()
    #>>> fastqcAnalysis(moveGzFiles = True)
    """

    cwd = os.getcwd()

    ## Fastq Analysis ##
    print("\n< FASTQC - Quality control analysis >")
    fastqgz_counter = 0
    for i in os.listdir(cwd):
        if 'fastq.gz' in i:
            fastqgz_counter += 1

    if fastqgz_counter == 0:
        print("\t> ERROR: No 'fastq.gz' files in current directory. Skipping fastqc analysis...")
    else:
        for i in os.listdir(cwd):
            if 'fastq.gz' in i:
                if seqPlatform == "MiSeq":
                    sample_idd, sep, rest = i.partition('_L001')
                    if sample_idd[-2:].isdigit() == True:
                        sample_true_id = sample_idd[:-4]
                    else:
                        sample_true_id = sample_idd[:-3]
                    print('\t> FastQC analysis running for sample: ', sample_true_id)
                    os.system('fastqc ' + i)
                elif seqPlatform == "NextSeq":
                    sample_true_id, sep, rest = i.partition('_R1_001')
                    print('\t> FastQC analysis running for sample: ', sample_true_id)
                    os.system('fastqc ' + i)

        print("\t> Fastqc Analysis...Done!")


def chrConcat(folder, deleteFaFiles=False, overwrite=False, mtChrAtEnd=True):
    """
    (path [, boolean, boolean]) -> full_karyo.fa

    Concatenates 'chrN.fa' files from folder into one 'full_karyo.fa' file. Path must be inside quotation marks.
    Files should be named like 'chr1.fa', 'chr2.fa', 'chrMT.fa', etc.
    Second and third arguments 'deleteFaFiles' / 'overwrite' are optional and are 'False' by default.

    Input:
    chr1.fa
    chr2.fa
    chrX.fa
    ...

    Output:
    full_karyo.fa

    Example:
    >>> chrConcat('/home/dfernandes/NGSdata/Indexes/chroms_hg38', deleteFaFiles=True)
    """
    print("\n< CHRCONCAT - Concatenating chromosomes >")

    ## Error checking
    if type(folder) != str:
        print("\t> ERROR: Path must be a 'str'. Concatenation aborted...")
        quit()
    else:
        pass

    try:
        os.chdir(folder)
    except FileNotFoundError:
        print("\t> ERROR: Specified path does not exist. Concatenation aborted...")
        quit()

    faCount = 0
    for i in os.listdir(folder):
        if i[-3:] == '.fa':
            faCount += 1

    if faCount == 0:
        print("\t> ERROR: No '.fa' files in the specified path. Concatenation aborted...")
        quit()
    else:
        pass

    if overwrite == False:
        for i in os.listdir(folder):
            if i == 'full_karyo.fa':
                print(
                    "\t> ERROR: Concatenated file 'full_karyo.fa' already exists. To enable overwriting, please use 'overwrite = True' when calling this function. Concatenation aborted...")
                quit()
            else:
                pass
    elif overwrite == True:
        pass

    ## If no errors occur, start here
    chr_list = []
    chr_index_list = []
    chr_index_list_ordered = []
    chr_index_list_double_digits = []
    os.chdir(folder)
    cwd = os.getcwd()
    file_list = os.listdir(cwd)

    ## Rejects files that have '_' in their names
    for i in file_list:
        if '_' in i:
            pass
        else:
            chr_list.append(i)

    ## Ordering chromosome files in numerical order
    for i in chr_list:
        index = i[3:]
        index = index.rstrip('.fa')
        chr_index_list.append(index)

    for i in chr_index_list:
        if len(i) == 1:
            chr_index_list_ordered.append(i)
        if len(i) == 2:
            chr_index_list_double_digits.append(i)

    chr_index_list_ordered.sort()
    chr_index_list_double_digits.sort()
    chr_index_list_ordered.extend(chr_index_list_double_digits)
    z = 0
    for i in chr_index_list_ordered:

        if i[0].isdigit() == False:
            pass
        else:
            chr_index_list_ordered.insert(z, chr_index_list_ordered.pop(chr_index_list_ordered.index(i)))
            z += 1

    ## Preparing str with final command
    chr_list.sort()
    command = ''
    for i in chr_index_list_ordered:
        for it in chr_list:
            if it[3:] == i + '.fa':
                command += it + ' '

    if mtChrAtEnd == True:
        chr_index_list_ordered += [chr_index_list_ordered.pop(chr_index_list_ordered.index("M"))]

    if command == '':
        print(
            "\t> ERROR: Fasta files naming is not correct. Please check 'help(chrConcat)' for more information. Concatenation aborted...")
        quit()
    else:
        print("\t> Concatenating chromosomes as:")
        print(chr_index_list_ordered)
        os.system('cat ' + command + ' > full_karyo.fa')

    if deleteFaFiles == True:
        cwd_list = os.listdir(cwd)
        for i in cwd_list:
            if 'full_karyo.fa' != i:
                os.remove(i)
    elif deleteFaFiles == False:
        pass

    print("\t> chromosomes' concatenation...Done!")


def makeIndexGenome(filePath, bwa=True, samtools=True, gatkpicard=True, overwrite=False):
    """
    (path [,boolean,boolean,boolean,boolean]) -> various_files

    Indexes the reference genome 'filePath' in FASTA format with a *.fa, *.fna, or *.fasta extension with BWA, SAMTOOLS, and GATK.
    Requires 'jarwrapper'.

    :param filePath:
    (str)   Path to the FASTA sequence to index

    :param bwa:
    (boolean)   Index with BWA. Default=True.

    :param samtools:
    (boolean)   Index with SAMTOOLS. Default=True.

    :param gatkpicard:
    (boolean)   Index with PICARD for GATK. Default=True.

    :param overwrite:
    (boolean)   Overwrite existing output files. Default=False.

    Example:
    >>> makeIndexGenome('/home/dfernandes/NGSdata/Indexes/chroms_hg38/full_karyo.fa', bwa=True, samtools=True, overwrite=True)
    """

    print("\n< BWA/SAMTOOLS - Indexing reference genome >")

    ## Error checking
    if type(filePath) != str:
        print("\t> ERROR: Path must be a 'str'. Indexing aborted...")
        quit()
    else:
        pass

    cwd, sep, id = filePath.rpartition('/')

    try:
        os.chdir(cwd)
    except FileNotFoundError:
        print("\t> ERROR: Specified folder does not exist. Indexing aborted...")
        quit()

    if os.path.isfile(filePath) == False:
        print("\t> ERROR: Specified file does not exist. Indexing aborted...")
        quit()

    os.chdir(cwd)
    if overwrite == False:
        if id + '.amb' in os.listdir(cwd) and id + '.ann' in os.listdir(cwd) and id + '.pac' in os.listdir(cwd):
            print(
                "\t> ERROR: Indexing files found in folder. To enable overwriting, please use 'overwrite = True' when calling this function. Indexing aborted...")
            quit()
    elif overwrite == True:
        pass

    ## If no errors occur, start here
    if bwa == True:
        print("\t> Indexing reference genome " + "'" + id + "'" + " with BWA in '" + cwd + "'...")
        print('\t\tbwa index -a bwtsw ' + filePath)
        os.system('bwa index -a bwtsw ' + filePath)
        print("\n",end='')

    if samtools == True:
        print("\t> Indexing reference genome " + "'" + id + "'" + " with SAMTOOLS in '" + cwd + "'...")
        print('\t\tsamtools faidx ' + filePath)
        os.system('samtools faidx ' + filePath)
        print("\n",end='')

    if gatkpicard == True:
        print("\t> Indexing reference genome " + "'" + id + "'" + " with PICARD for GATK in '" + cwd + "'...")
        patho, sep, extensiono = filePath.rpartition('.')
        patho2, sep2, extensiono2 = patho.rpartition('/')
        if overwrite == True and extensiono2 + '.dict' in os.listdir(cwd):
            os.remove(extensiono2 + '.dict')
        print('\t\tpicard.jar CreateSequenceDictionary R= ' + filePath + ' O= ' + patho + '.dict')
        os.system('picard.jar CreateSequenceDictionary R= ' + filePath + ' O= ' + patho + '.dict')
        print("\n",end='')
    print("\t> Reference genome indexing... Done!")


def cutadapt(nThreads = 1):
    """
    Removes the adapter sequences from '.fastq.gz' or '.fastq' files in current folder.
    Tested and working for version 1.5. More recent versions such as 1.8 include commas as thousand separators, which breaks the "End_analysis" function.
    IMPORTANT: File names cannot include the character %.

    Input:
    bar_S10_L001_R1_001.fastq.gz
    OR
    bar_S10_L001_R1_001.fastq

    Output:
    bar%S10_L001_R1_001.trimmed.fastq
    bar%S10_L001_R1_001.trimmed.log

    Example:
    >>> cutadapt()
    """
    cwd = os.getcwd()
    # print("> Current working directory is: ",cwd)
    print("\n< CUTADAPT - Removing adapters >")
    fout.write("\n< CUTADAPT - Removing adapters >\n")
    fastqGzCount = 0
    fastqCount = 0
    nThreads = str(nThreads)
    print(
        "ATTENTION: Cutadapt for paired-end data, as of now, requires the user to manually edit the names of the input 'fastq.gz' file as 'NAME_L1_P1.fastq.gz', 'NAME_L1_P2.fastq.gz', etc.")
    # print("ATTENTION: Cutadapt for paired-end data, as of now, only works for one sample per folder. Please create a folder for each sample.")

    if seqType == "singleEnd":
        for i in os.listdir(cwd):
            if 'fastq.gz' in i:
                fastqGzCount += 1
                if seqPlatform == "MiSeq":
                    sample_idd, sep, rest = i.partition('_L001')
                    if sample_idd[-2:].isdigit() == True:
                        sample_id = sample_idd[:-4]
                    else:
                        sample_id = sample_idd[:-3]
                elif seqPlatform == "NextSeq":
                    sample_id, sep, rest = i.partition('_R1_001')

                ## For MiSeq
                ## Deletes duplicate normal '*.fastq' file if one with the divider '%' already exists
                if seqPlatform == "MiSeq":
                    ############ This first statement might be bugged
                    print(sample_id + '%' + rest)
                    directory = os.listdir(cwd)
                    print("Renaming", i)
                    os.rename(directory[directory.index(i)], sample_id + '%' + rest)

                ## For NextSeq
                ## Deletes duplicate normal '*.fastq' file if one with the divider '%' already exists
                elif seqPlatform == "NextSeq":
                    print("I=",i)
                    print(sample_id + '%' + sep + rest)
                    if "%" not in i:
                        print("HERE")
                        print("Renaming", i)
                        directory = os.listdir(cwd)
                        os.rename(directory[directory.index(i)], sample_id + '%' + sep + rest)

        for i in os.listdir(cwd):
            # if '.fastq' in i and '.trimmed' not in i:
            if '.fastq.gz' in i and '.trimmed' not in i and "%" in i:
                # fastqCount += 1
                s, sep, r = i.partition('%')
                ## Adapter sequence is the sequence of the IS3 adapter (13bp), extended with the remaining complement from IS2, totalling 33bp
                ## As a note, IS1 and IS2 are inverted complements
                ## IS3     5? AGATCGGAAGAGC ACACGTCTGAACTCCAGTCAC 3?
                ## IS2     3? TCTAGCCTTCTCG TGTGCAGACTTGAGGTCAGTG 5?
                if sep == '%':
                    fastqCount += 1
                    print("\t> Cutting adapters from file " + "'" + i + "'" + "...")
                    fout.write("\t> Cutting adapters from file " + "'" + i + "'" + "...\n")
                    print('\tcutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -j ' + nThreads +' -m 17 ' + i + ' > ' + i[:-9] + '.trimmed.log' + "\n")
                    fout.write(time.strftime("\t%H:%M:%S %Z") + '\tcutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -j ' + nThreads + ' -m 17 ' + i + ' > ' + i[:-9] + '.trimmed.fastq' + ' 2> ' + i[:-9] + '.trimmed.log' + "\n\n")
                    os.system('cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 1 -j ' + nThreads + ' -m 17 ' + i + ' > ' + i[:-9] + '.trimmed.fastq' + ' 2> ' + i[:-9] + '.trimmed.log')

#    ISSUES HERE!!! UP


                    # fout = open(s+'_commands.txt', 'w') #s makes one file per sample
                    # with open(sample_id+'_commands.txt', 'w') as fout:

    elif seqType == "pairedEnd":
        listSamps = []
        for i in os.listdir(cwd):
            if '.fastq.gz' in i and '.trimmed' not in i:
                s, r = re.split("_P\d", i)
                listSamps.append(s)
        for i in set(listSamps):
            #            if '.fastq.gz' in i and '.trimmed' not in i:
            # fastqCount += 1
            # print(i)
            s1 = str(i + "_P1.fastq.gz")
            s2 = str(i + "_P2.fastq.gz")
            # s, sep, r = i.partition('_P1')
            # s, sep, r = re.split("_P1\d",i,1)
            ## Adapter sequence is the sequence of the IS3 adapter (13bp), extended with the remaining complement from IS2, totalling 33bp
            ## As a note, IS1 and IS2 are inverted complements
            ## IS3     5? AGATCGGAAGAGC ACACGTCTGAACTCCAGTCAC 3?
            ## IS2     3? TCTAGCCTTCTCG TGTGCAGACTTGAGGTCAGTG 5?
            #                if sep == '%':
            if s1 and s2 in os.listdir(cwd):
                fastqCount += 1
                print("\t> Cutting adapters from file " + "'" + i + "'" + "...")
                fout.write("\t> Cutting adapters from file " + "'" + i + "'" + "...\n")
                print(
                    '\tcutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGC -O 1 -m 17 -o ' + i + "P1%trimmed.fastq -p " + i + "P2%trimmed.fastq " + i + "_P1.fastq.gz " + i + "_P1.fastq.gz > " + i + ".trimmed.log" + "\n")
                fout.write(time.strftime(
                    "\t%H:%M:%S %Z") + '\tcutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGC -O 1 -m 17 -o ' + i + "P1%trimmed.fastq -p " + i + "P2%trimmed.fastq " + i + "_P1.fastq.gz " + i + "_P1.fastq.gz > " + i + ".trimmed.log" + "\n\n")
                os.system(
                    'cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGC -O 1 -m 17 -o ' + i + "P1%trimmed.fastq -p " + i + "P2%trimmed.fastq " + i + "_P1.fastq.gz " + i + "_P1.fastq.gz > " + i + ".trimmed.log")

                # fout = open(s+'_commands.txt', 'w') #s makes one file per sample
                # with open(sample_id+'_commands.txt', 'w') as fout:

    # os.chdir(originalWD)
    # cwd = os.getcwd()
    # print("\n\t< Collecting data... >")
    # fout.write("\n\t< Collecting data... >")
    if seqType == "pairedEnd":
        print("\n\t< 'Data collection' not available for paired end reads...Yet. >")
        fout.write("\n\t< 'Data collection' not available for paired end reads...Yet. >")

    elif seqType == "singleEnd":
        cutout = open('Cutadapt_Summary.txt', 'w')
        cutout.write("Sample\tTotal_reads\tTrimmed_reads\tNum_lanes\n")

        ## For MiSeq
        if seqPlatform == "MiSeq":
            sampleDict = {}
            for i in os.listdir(cwd):
                #print(i)
                # ja = re.search(r'_S\d\d_L\d\d\d\.q30\.rmdup\.bam$', i)
                ja = re.search(r'%_R1_\d\d\d\.trimmed\.log$', i)
                ba = re.search(r'_S\d{1,2}_L\d\d\d%_R1_\d\d\d\.', i)
                if ja is not None and ba is None:
                    sampleName = i[:ja.start()]
                    #print(sampleName)

                    if sampleName in sampleDict.keys():
                        if i in sampleDict[sampleName]:
                            pass

                        else:
                            sampleDict[sampleName].append(i)
                            #print("b")
                    else:
                        sampleDict[sampleName] = [i]
        ## For NextSeq
        elif seqPlatform == "NextSeq":
            sampleDict = {}
            for i in os.listdir(cwd):
                #print(i)
                # ja = re.search(r'_S\d\d_L\d\d\d\.q30\.rmdup\.bam$', i)
                ja = re.search(r'_S\d{1,2}_L\d\d\d%_R1_\d\d\d\.trimmed\.log$', i)
                #print(ja)
                #ba = re.search(r'\.trimmed\.log$', i)
                #if ja is not None and ba is None:
                if ja is not None:
                    sampleName = i[:ja.start()]
                    print(sampleName)

                    if sampleName in sampleDict.keys():
                        if i in sampleDict[sampleName]:
                            pass

                        else:
                            sampleDict[sampleName].append(i)
                            print("b")
                    else:
                        sampleDict[sampleName] = [i]

            #print(sampleDict)

        for i in sampleDict.keys():
            print("\n\t< Collecting data for sample "+i+" >")
            fout.write("\n\t< Collecting data for sample "+i+" >")
            procReadsTot = []
            trimReadsTot = []
            for a in sampleDict[i]:
                print("\tFrom " + a + "...")
                fout.write("\tFrom  " + a + "...")
                ## Processed and trimmed reads
                with open(a, 'r') as fin:
                    for line in fin:
                        ## Processed reads
                        if "Total reads processed:" in line:
                            text = line.split()
                            ba = text[3]
                            ba = ba.replace(',', '')
                            procReadsTot.append(int(ba))
                        ## Trimmed reads
                        if "Reads with adapters:" in line:
                            text = line.split()
                            ca = text[3]
                            ca = ca.replace(',', '')
                            trimReadsTot.append(int(ca))
                            #print("Trimmed reads: "+ca)
            cutout.write(i + "\t" + str(sum(procReadsTot)) + "\t" + str(sum(trimReadsTot)) + "\t" + str(len(procReadsTot)) + "\n")


    if fastqGzCount == 0 and fastqCount == 0:
        print("\n\t> ERROR: No '.fastq.gz' or '.fastq' files in folder. Cutadapt aborted...")
        fout.write(time.strftime(
            "\n\t%H:%M:%S %Z ") + "> ERROR: No '.fastq.gz' or '.fastq' files in folder. Cutadapt aborted...\n")
    else:
        print("\t> Adatper sequences removal... Done!\n\t-----")
        fout.write("\t> Adatper sequences removal... Done!\n\t-----\n")


def trimFastq(bp):
    # def trimFastq(group_tag):
    """
    Removes the 2 base pairs at both ends of each read.
    WARNING: Currently breaks the downstream analysis as output is not automatically recognized.

    :param bp:
    (int) number of bases to trim at each end of each read.

    Input:
    foo%S10_L001_R1_001.trimmed.fastq

    Output:
    foo.seqtk.fastq
    foo.seqtk.log

    Example:
    >>> trimFastq(2)
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    trimmedCount = 0

    for i in os.listdir(cwd):
        if re.search("trimmed.fastq", i) != None:
            trimmedCount += 1
            sample_id, rest = i.split('%')
            print("\t> Trimming terminal bases from sequences from " + "'" + i + "'" + "...")
            # print('samtools bam2fq '+i+' > '+sample_id+'.fastq')
            # os.system('samtools bam2fq '+i+' > '+sample_id+'.fastq')
            print('\tseqtk trimfq -b ' + str(bp) + ' -e ' + str(
                bp) + ' ' + i + ' > ' + sample_id + '%seqtk.fastq' + ' 2> ' + sample_id + '%seqtk.log' + "\n")
            os.system('seqtk trimfq -b ' + str(bp) + ' -e ' + str(
                bp) + ' ' + i + ' > ' + sample_id + '%seqtk.fastq' + ' 2> ' + sample_id + '%seqtk.log')
            fout.write(
                "\n\t> Trimming terminal bases from sequences from " + "'" + sample_id + "'" + "..." + time.strftime(
                    "\n\t%H:%M:%S %Z ") + '\tseqtk trimfq -b ' + str(bp) + ' -e ' + str(
                    bp) + ' ' + i + ' > ' + sample_id + '%seqtk.fastq' + ' 2> ' + sample_id + '%seqtk.log' + "\n")
            seqtkSim = 1

    if trimmedCount == 0:
        print("\t> ERROR: No 'trimmed.fastq' files in folder. No sequences to trim. Trimming aborted...")
        fout.write(time.strftime(
            "\n\t%H:%M:%S %Z ") + "> ERROR: No 'trimmed.fastq' files in folder. No sequences to trim. Trimming aborted...\n")

    else:
        print("\t> Terminal bases trimming... Done!\n\t-----")
        fout.write("\n\t> Terminal bases trimming... Done!\n\t-----\n")


def bwaAlign(ref_gen, nThreads = 1, overwrite=False):
    """
    (path [,int,boolean]) -> various_files

    Aligns the 'trimmed.fastq' files in current folder to the indexed genome 'ref_gen'. Path must be inside quotation marks.
    Argument 'overwrite' is optional and 'False' by default.

    :param ref_gen:
    (str) path to reference genome.

    :param nThreads:
    (int) number of CPU execution threads.

    :param overwrite:
    (boolean) overwrite existing *.sai files?

    Input:
    foo%S10_L001_R1_001.trimmed.fastq

    Output:
    foo.sai
    foo.bwa_aln.log

    Example:
    >>> bwaAlign(ref_gen='/home/dfernandes/NGSdata/Indexes/chroms_hg38/full_karyo.fa', nThreads=6, overwrite=True)
    """

    cwd = os.getcwd()
    # print("> Current working directory is: ",cwd)
    print("\n< BWA - Aligning to ref genome >")
    fout.write("\n< BWA - Aligning to ref genome >")

    if seqtkSim == 1:
        fileExp = "seqtk.fastq$"
    else:
        fileExp = "trimmed.fastq$"

    ## Aligning the cutadapt-trimmed files to the reference genome
    trimmedCount = 0
    nThreads = str(nThreads)
    ##### TRY NEW IDENTIFIER!
    for i in os.listdir(cwd):
        if re.search(fileExp, i) != None:
            trimmedCount += 1
            sample_id, rest = i.split('%')

            ## Aligning the cutadapt-trimmed files to the reference genome
            ## Checks if there is already a '.sai' file for the correspondent iterator and 'overwrites' it or not
            if sample_id + '.sai' in os.listdir(cwd):
                if overwrite == False:
                    pass
                elif overwrite == True:
                    print("\t> Aligning trimmed sequences from " + "'" + i + "'" + " to reference genome...")
                    print(
                        '\tbwa aln -l 1000 -t ' + nThreads + ' ' + ref_gen + ' ' + i + ' > ' + sample_id + '.sai' + ' 2> ' + sample_id + '.bwa_aln.log' + "\n")
                    os.system(
                        'bwa aln -l 1000 -t ' + nThreads + ' ' + ref_gen + ' ' + i + ' > ' + sample_id + '.sai' + ' 2> ' + sample_id + '.bwa_aln.log')
                    fout.write(
                        "\n\t> Aligning trimmed sequences from " + "'" + i + "'" + " to reference genome..." + '\n\tbwa aln -l 1000 -t ' + nThreads + ' ' + ref_gen + ' ' + i + ' > ' + sample_id + '.sai' + ' 2> ' + sample_id + '.bwa_aln.log' + "\n")

            ## If there is no '.sai' file for this iterator, run normally
            else:
                print("\t> Aligning trimmed sequences from " + "'" + i + "'" + " to reference genome...")
                print(
                    '\tbwa aln -l 1000 -t ' + nThreads + ' ' + ref_gen + ' ' + i + ' > ' + sample_id + '.sai' + ' 2> ' + sample_id + '.bwa_aln.log' + "\n")
                os.system(
                    'bwa aln -l 1000 -t ' + nThreads + ' ' + ref_gen + ' ' + i + ' > ' + sample_id + '.sai' + ' 2> ' + sample_id + '.bwa_aln.log')
                fout.write(
                    "\n\t> Aligning trimmed sequences from " + "'" + i + "'" + " to reference genome..." + time.strftime(
                        "\n\t%H:%M:%S %Z ") + '\tbwa aln -l 1000 -t ' + nThreads + ' ' + ref_gen + ' ' + i + ' > ' + sample_id + '.sai' + ' 2> ' + sample_id + '.bwa_aln.log' + "\n")

            if sample_id + '.bwa_aln.log' in os.listdir(cwd):
                with open(sample_id + '.bwa_aln.log', 'r') as fin_log:
                    file_log = fin_log.readlines()
                    if any("fail to locate the index" in s for s in file_log):
                        print(
                            "\t> ERROR: Failed to locate specified index file " + ref_gen + ". Aborting downstream analysis.")
                        fout.write(time.strftime(
                            "\n\t%H:%M:%S %Z ") + "\t> ERROR: Failed to locate specified index file " + ref_gen + ". Aborting downstream analysis.")
                        exit()
            else:
                print("\t> Aligning to reference genome... Done!\n\t-----")
                fout.write("\n\t> Aligning to reference genome... Done!\n\t-----\n")

    if trimmedCount == 0:
        print("\t> ERROR: No 'trimmed.fastq' files in folder. No sequences to align. Aligning aborted...")
        fout.write(time.strftime(
            "\n\t%H:%M:%S %Z ") + "> ERROR: No 'trimmed.fastq' files in folder. No sequences to align. Aligning aborted...\n")


def tagAndConvertSb(ref_gen, group_tag = "Unnamed", qualThreshold = 30):
    """
    (path [,str,int]) -> various_files

    Adds read groups tags and converts SAM to BAM.

    :param ref_gen:
    (str) path to reference genome.

    :param group_tag:
    (str) read groups tag for @RG\tID.

    :param qualThreshold:
    (int) quality score (QC) filter.

    Input:
    foo%S10_L001_R1_001.trimmed.fastq
    foo.sai

    Output:
    foo.q30.sam
    foo.q30.bam

    Example:
    ########>>> bwaAlign('/home/dfernandes/NGSdata/Indexes/chroms_hg38/full_karyo.fa', overwrite=True)
    """

    ### insert 'overwrite' option and file count ERROR


    cwd = os.getcwd()
    print("\n< BWA - Tagging and converting SAM>BAM >")
    fout.write("\n< BWA - Tagging and converting SAM>BAM >")

    if seqtkSim == 1:
        fileExp = "seqtk.fastq$"
    else:
        fileExp = "trimmed.fastq$"
    qualThreshold = str(qualThreshold)


    if seqType == "singleEnd":
        ## Converting '*.sai' to '*.sam' and sequently '*.bam', and adding tags
        for i in os.listdir(cwd):
            if re.search(".sai$", i) != None:
                sample_id = i[:-4]

                for it in os.listdir(cwd):
                    if '%' in it and sample_id in it and re.search(fileExp, it) != None:
                        file_id, sep, rest = it.rpartition('%')
                        if seqtkSim == 0:
                            if file_id == sample_id and file_id in it and 'trimmed.fastq' in it:
                                print("\t> Converting " + i + " to '.sam' and adding tags...")
                                print(
                                    "\t" + r"bwa samse -r '@RG\tID:" + group_tag + r"\tSM:" + sample_id + r"\tPL:ILLUMINA' " + ref_gen + ' ' + i + ' ' + it + ' > ' + sample_id + '.sam')
                                os.system(
                                    r"bwa samse -r '@RG\tID:" + group_tag + r"\tSM:" + sample_id + r"\tPL:ILLUMINA' " + ref_gen + ' ' + i + ' ' + it + ' > ' + sample_id + '.sam')
                                fout.write("\n\t> Converting " + i + " to '.sam' and adding tags..." + time.strftime(
                                    "\n\t%H:%M:%S %Z ") + "\t" + r"bwa samse -r '@RG\tID:" + group_tag + r"\tSM:" + sample_id + r"\tPL:ILLUMINA' " + ref_gen + ' ' + i + ' ' + it + ' > ' + sample_id + '.sam' + '\n')
                                print("\t> Converting " + sample_id + '.sam' + " to '.q" + qualThreshold + ".bam'...")
                                print(
                                    "\t"r"samtools view -Sb -q "+ qualThreshold + " -F 4 " + sample_id + ".sam" + " > " + sample_id + '.q'+ qualThreshold + '.bam')
                                os.system(
                                    r"samtools view -Sb -q "+ qualThreshold + " -F 4 " + sample_id + ".sam" + " > " + sample_id + '.q'+ qualThreshold + '.bam')
                                fout.write("\t> Converting " + sample_id + '.sam' + " to '.q" + qualThreshold + ".bam'..." + time.strftime(
                                    "\n\t%H:%M:%S %Z ") + "\t"r"samtools view -Sb -q "+ qualThreshold + " -F 4 " + sample_id + ".sam" + " > " + sample_id + '.q'+ qualThreshold + '.bam' + '\n')
                                print("\t")
                        if seqtkSim == 1:
                            if file_id == sample_id and file_id in it and 'seqtk.fastq' in it:
                                print("\t> Converting " + i + " to '.sam' and adding tags...")
                                print(
                                    "\t" + r"bwa samse -r '@RG\tID:" + group_tag + r"\tSM:" + sample_id + r"\tPL:ILLUMINA' " + ref_gen + ' ' + i + ' ' + it + ' > ' + sample_id + '.sam')
                                os.system(
                                    r"bwa samse -r '@RG\tID:" + group_tag + r"\tSM:" + sample_id + r"\tPL:ILLUMINA' " + ref_gen + ' ' + i + ' ' + it + ' > ' + sample_id + '.sam')
                                fout.write("\n\t> Converting " + i + " to '.sam' and adding tags..." + time.strftime(
                                    "\n\t%H:%M:%S %Z ") + "\t" + r"bwa samse -r '@RG\tID:" + group_tag + r"\tSM:" + sample_id + r"\tPL:ILLUMINA' " + ref_gen + ' ' + i + ' ' + it + ' > ' + sample_id + '.sam' + '\n')
                                print("\t> Converting " + sample_id + '.sam' + " to '.q" + qualThreshold + ".bam'...")
                                print(
                                    "\t"r"samtools view -Sb -q "+ qualThreshold + " -F 4 " + sample_id + ".sam" + " > " + sample_id + '.q'+ qualThreshold + '.bam')
                                os.system(
                                    r"samtools view -Sb -q "+ qualThreshold + " -F 4 " + sample_id + ".sam" + " > " + sample_id + '.q'+ qualThreshold + '.bam')
                                fout.write("\t> Converting " + sample_id + '.sam' + " to '.q" + qualThreshold + ".bam'..." + time.strftime(
                                    "\n\t%H:%M:%S %Z ") + "\t"r"samtools view -Sb -q "+ qualThreshold + " -F 4 " + sample_id + ".sam" + " > " + sample_id + '.q'+ qualThreshold + '.bam' + '\n')
                                print("\t")


    elif seqType == "pairedEnd":
        listSamps = []
        for i in os.listdir(cwd):
            if i[-4:] == ".sai":
                s, r = re.split("_P\d", i)
                listSamps.append(s)
        for i in set(listSamps):
            s1 = str(i + "_P1.sai")
            s2 = str(i + "_P2.sai")
            if seqtkSim == 0:
                t1 = str(i + "_P1%trimmed.fastq")
                t2 = str(i + "_P2%trimmed.fastq")
            elif seqtkSim == 1:
                t1 = str(i + "_P1%seqtk.fastq")
                t2 = str(i + "_P2%seqtk.fastq")
            if s1 and s2 and t1 and t2 in os.listdir(cwd):
                print("\t> Converting " + i + " to '.sam' and adding tags...")
                print(
                    "\t" + r"bwa sampe -r '@RG\tID:" + group_tag + r"\tSM:" + i + r"\tPL:ILLUMINA' " + ref_gen + ' ' + s1 + ' ' + s2 + ' ' + t1 + ' ' + t2 + ' > ' + i + '.sam')
                os.system(
                    r"bwa sampe -r '@RG\tID:" + group_tag + r"\tSM:" + i + r"\tPL:ILLUMINA' " + ref_gen + ' ' + s1 + ' ' + s2 + ' ' + t1 + ' ' + t2 + ' > ' + i + '.sam')
                fout.write("\n\t> Converting " + i + " to '.sam' and adding tags..." + time.strftime(
                    "\n\t%H:%M:%S %Z ") + "\t" + r"bwa sampe -r '@RG\tID:" + group_tag + r"\tSM:" + i + r"\tPL:ILLUMINA' " + ref_gen + ' ' + s1 + ' ' + s2 + ' ' + t1 + ' ' + t2 + ' > ' + i + '.sam' + '\n')
                print("\t> Converting " + i + '.sam' + " to '.q" + qualThreshold + ".bam'...")
                print("\t"r"samtools view -Sb -q "+ qualThreshold + " -F 4 " + i + ".sam" + " > " + i + '.q'+ qualThreshold + '.bam')
                os.system(r"samtools view -Sb -q "+ qualThreshold + " -F 4 " + i + ".sam" + " > " + i + '.q'+ qualThreshold + '.bam')
                fout.write("\t> Converting " + i + '.sam' + " to '.q" + qualThreshold + ".bam'..." + time.strftime(
                    "\n\t%H:%M:%S %Z ") + "\t"r"samtools view -Sb -q "+ qualThreshold + " -F 4 " + i + ".sam" + " > " + i + '.q'+ qualThreshold + '.bam' + '\n')
                print("\t")

    print("\t> Tagging and converting... Done!")
    fout.write("\n\t> Tagging and converting... Done!\n\t-----\n")


def samToolsSort(version="1.3.1"):
    """
    ([str]) -> various_files

    Sorts the BAM files.

    :param version:
    (str) samtools' version.

    Input:
    foo.q30.bam

    Output:
    foo.q30.sort.bam

    Example:
    >>> bwaAlign('/home/dfernandes/NGSdata/Indexes/chroms_hg38/full_karyo.fa', overwrite=True)
    """

    cwd = os.getcwd()
    print("\n< SAMTOOLS - Sorting BAM files >")
    fout.write("\n< SAMTOOLS - Sorting BAM files >")

    if version[0:3] == "0.1":
        for i in os.listdir(cwd):
            if '.bam' == i[-4:] and '.sort' not in i and '.rmdup' not in i:
                sample_id = i[:-4]
                print("\t> Sorting " + i + "...")
                print("\t" + 'samtools sort ' + i + ' ' + sample_id + '.sort' + "\n")
                os.system('samtools sort ' + i + ' ' + sample_id + '.sort')
                fout.write("\n\t> Sorting " + i + "..." + time.strftime("\n\t%H:%M:%S %Z ") + '\tsamtools sort ' + i + ' ' + sample_id + '.sort' + "\n")

    if version == "1.3.1":
        for i in os.listdir(cwd):
            if '.bam' == i[-4:] and '.sort' not in i and '.rmdup' not in i:
                sample_id = i[:-4]
                print("\t> Sorting " + i + "...")
                print("\t" + 'samtools sort ' + i + ' -o ' + sample_id + '.sort.bam' + "\n")
                os.system('samtools sort ' + i + ' -o ' + sample_id + '.sort.bam')
                fout.write("\n\t> Sorting " + i + "..." + time.strftime("\n\t%H:%M:%S %Z ") + '\tsamtools sort ' + i + ' -o ' + sample_id + '.sort.bam' + "\n")


    print("\t> Sorting BAM files... Done!")
    fout.write("\n\t> Sorting BAM files... Done!\n\t-----\n")


def fastqcBeforeRmdup():
    """
    () -> html,zip

    For current folder, executes FastQC analysis for each 'fastq.gz' file, creating individual folders.

    Input:
    foo.q30.sort.bam

    Output:
    foo.q30.sort.zip
    foo.q30.sort.html

    Example:
    >>> fastqcBeforeRmdup()
    """

    cwd = os.getcwd()

    ## Fastq Analysis ##
    print("\n< FASTQC - Quality control analysis >")
    fastq_counter = 0
    for i in os.listdir(cwd):
        if '.sort.bam' in i:
            fastq_counter += 1

    if fastq_counter == 0:
        print("\t> ERROR: No '.sort.bam' files in current directory. Skipping fastqc analysis...")
    else:
        for i in os.listdir(cwd):
            if '.sort.bam' == i[-9:]:
                sample_id = i[:-9]
                print('\t> FastQC analysis running for sample: ', sample_id)
                os.system('fastqc ' + i)


        print("\t> Fastqc Analysis...Done!")


def rmvDups(picard=False):
    """
    ([boolean]) -> various_files

    Removes duplicated sequences using SAMTOOLS or PICARD/GATK.

    Input:
    foo.q30.sort.bam

    Output:
    foo.q30.rmdup.bam

    Example:
    ######>>> bwaAlign('/home/dfernandes/NGSdata/Indexes/chroms_hg38/full_karyo.fa', overwrite=True)
    """

    cwd = os.getcwd()
    print("\n< SAMTOOLS - Removing duplicates >")
    fout.write("\n< SAMTOOLS - Removing duplicates >")

    ## Selects whether to use Picard or Samtools
    if picard == True:
        for i in os.listdir(cwd):
            if '.sort.bam' == i[-9:]:
                sample_id = i[:-9]
                print("\t> Removing duplicates using Picard from file " + i + "...")
                print(
                    "\t" + 'picard.jar MarkDuplicates I=' + i + ' O=' + sample_id + '.rmdup.bam METRICS_FILE=' + sample_id + '_metrics.txt')
                os.system(
                    'picard.jar MarkDuplicates I=' + i + ' O=' + sample_id + '.rmdup.bam METRICS_FILE=' + sample_id + '_metrics.txt')
                fout.write("\n\t> Removing duplicates from file " + i + "using PicardTools..." + time.strftime(
                    "\n\t%H:%M:%S %Z ") + "\t" + 'picard.jar MarkDuplicates I=' + i + ' O=' + sample_id + '.rmdup.bam METRICS_FILE=' + sample_id + '_metrics.txt' + '\n')
                print("\t")

    elif picard == False:
        for i in os.listdir(cwd):
            if '.sort.bam' == i[-9:]:
                sample_id = i[:-9]
                print("\t> Removing duplicates from file " + i + "...")
                print("\t" + 'samtools rmdup -s ' + i + ' ' + sample_id + '.rmdup.bam')
                os.system('samtools rmdup -s ' + i + ' ' + sample_id + '.rmdup.bam')
                fout.write("\n\t> Removing duplicates from file " + i + "..." + time.strftime(
                    "\n\t%H:%M:%S %Z ") + "\t" + 'samtools rmdup -s ' + i + ' ' + sample_id + '.rmdup.bam' + '\n')
                print("\t")

    for i in os.listdir(cwd):
        if '.rmdup.bam' == i[-10:]:
            sample_id = i[:-10]
            os.system('samtools flagstat ' + i + ' > ' + sample_id + '.summary.txt')
            print("\t> Generating 'flagstats' summary file for " + i + "...")
            print('\tsamtools flagstat ' + i + ' > ' + sample_id + '.summary.txt\n')
            fout.write("\n\t> Generating 'flagstats' summary file for " + i + "..." + time.strftime(
                "\n\t%H:%M:%S %Z ") + '\tsamtools flagstat ' + i + ' > ' + sample_id + '.summary.txt\n')

    print("\t> Removing duplicates... Done!")
    fout.write(("\n\t> Removing duplicates... Done!\n\t-----\n"))


def makeIndexBam():
    """

    Input:
    foo.q30.rmdup.bam

    Output:
    foo.q30.rmdup.bam.bai
    """

    cwd = os.getcwd()
    print("\n< SAMTOOLS - Indexing BAM files >")
    fout.write(("\n< SAMTOOLS - Indexing BAM files >"))

    for i in os.listdir(cwd):
        if '.rmdup.bam' == i[-10:]:
            sample_id = i[:-10]
            print("\t> Indexing '*.bam' file {}...".format(sample_id))
            os.system('samtools index ' + i + ' ' + i + '.bai')
            fout.write("\n\t> Indexing '*.bam' file {}...".format(sample_id))
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + '\tsamtools index ' + i + ' ' + i + '.bai\n')

    print("\n\t> BAM file indexing... Done!")
    fout.write("\n\t> BAM file indexing... Done!\n\t-----\n")


def mapdamage(filepath, stats=False):
    """
    Runs the mapDamage code to track and quantify damage patterns in ancient DNA sequences.

    Input:
    foo.q30.rmdup.bam

    Output:
    foo (folder)
    Example:
    >>> mapdamage('/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa')

    """
    cwd = os.getcwd()
    print("\n< MAPDAMAGE - Damage patterns in ancient DNA >")
    fout.write("\n< MAPDAMAGE - Damage patterns in ancient DNA >")
    q30rmdup_counter = 0
    for i in os.listdir(cwd):
        if '.rmdup.bam' in i:
            q30rmdup_counter += 1

    if q30rmdup_counter == 0:
        print("\t> ERROR: No '.rmdup.bam' files in current directory. Aborting MAPDAMAGE analysis...")
        fout.write(time.strftime(
            "\n\t%H:%M:%S %Z ") + "> ERROR: No '.rmdup.bam' files in current directory. Aborting MAPDAMAGE analysis...")
        return None
    else:
        for i in os.listdir(cwd):
            if '.rmdup.bam' == i[-10:]:
                sample_id = i[:-14]
                print("\t> Performing MAPDAMAGE analysis for sample " + i + "...")
                if stats == False:
                    print(
                        "\t""mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id + " --no-stats" + "\n")
                    os.system(
                        "mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id + " --no-stats")
                    fout.write("\t> Performing MAPDAMAGE analysis for sample " + i + "..." + time.strftime(
                        "\n\t%H:%M:%S %Z ") + "\t""mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id + " --no-stats" + "\n")
                elif stats == True:
                    print(
                        "\t""mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id + "\n")
                    os.system(
                        "mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id)
                    fout.write("\t> Performing MAPDAMAGE analysis for sample " + i + "..." + time.strftime(
                        "\n\t%H:%M:%S %Z ") + "\t""mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id + "\n")
    print("\n\t> Damage patterns analysis... Done!")
    fout.write("\n\n\t> Damange patterns analysis... Done!\n\t-----\n")


def rescaleMapDamage(filepath):
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< MAPDAMAGE RESCALE - Rescale base quality scores with mapDamage algo >")
    fout.write("\n< MAPDAMAGE RESCALE - Rescale base quality scores with mapDamage algo >")

    # if seqtkSim == 1:
    #     fileExp = "seqtk.bam$"
    # else:
    #     fileExp = "q30.rmdup.bam$"

    for i in os.listdir(cwd):
        if re.search("q..\.rmdup\.bam$", i) != None:
            sample_id, sep, rest = i.partition('.bam')
            sample_id_true, rest = re.split('\.q..',i)
            sep = re.findall('\.q..',i)[0]
            if os.path.exists(cwd + "/mapdamage/mapdamage_" + sample_id_true) == True:
                if 'Stats_out_MCMC_correct_prob.csv' in os.listdir(cwd + "/mapdamage/mapdamage_" + sample_id_true):
                    print("mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id_true + " --rescale-only")
                    os.system("mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id_true + " --rescale-only")
                elif 'Stats_out_MCMC_correct_prob.csv' not in os.listdir(cwd + "/mapdamage/mapdamage_" + sample_id_true):
                    print("mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id_true + " --stats-only")
                    os.system("mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id_true + " --stats-only")
                    print("mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id_true + " --rescale-only")
                    os.system("mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id_true + " --rescale-only")
            else:
                print("mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id_true + " --rescale")
                os.system("mapDamage -i " + i + " -r " + filepath + " -d " + cwd + "/mapdamage/mapdamage_" + sample_id_true + " --rescale")
            mapDamageRescaleSim = 1


def fastqcAfterRmdup():
    """
    () -> html,zip

    For current folder, executes FastQC analysis for each 'fastq.gz' file, creating individual folders.

    Input:
    foo.q30.rmdup.bam

    Output:
    foo.q30.rmdup.zip
    foo.q30.rmdup.html

    Example:
    >>> fastqcAfterRmdup()
    """

    cwd = os.getcwd()

    ## Fastq Analysis ##
    print("\n< FASTQC - Quality control analysis >")
    fastq_counter = 0
    for i in os.listdir(cwd):
        if '.rmdup.bam' in i:
            fastq_counter += 1

    if fastq_counter == 0:
        print("\t> ERROR: No '.rmdup.bam' files in current directory. Skipping fastqc analysis...")
    else:
        for i in os.listdir(cwd):
            if '.rmdup.bam' == i[-10:]:
                sample_id = i[:-10]
                print('\t> FastQC analysis running for sample: ', sample_id)
                os.system('fastqc ' + i)

        print("\t> Fastqc Analysis...Done!")


def idxstats():
    """
    () -> txt

    Outputs a file with the distribution of the reads per chromosome for all files with pattern in folder.

    Input:
    bar.q30.rmdup.bam

    Output:
    bar.chr.txt

    Example:
    >>> idxstats()
    """
    cwd = os.getcwd()
    print("\n< IDXSTATS - Distribution of reads per chromosome >")
    fout.write("\n< IDXSTATS - Distribution of reads per chromosome >")

    for i in os.listdir(cwd):
        if '.rmdup.bam' == i[-10:]:
            sample_id = i[:-14]
            print("\t> Reading distribution of reads per chromosome from '" + i + "'...")
            print("\t" + r"samtools idxstats " + i + " > " + sample_id + ".chr.txt\n")
            os.system(r"samtools idxstats " + i + " > " + sample_id + ".chr.txt")
            fout.write("\n\t> Reading distribution of reads per chromosome from '" + i + "'..." + time.strftime(
                "\n\t%H:%M:%S %Z ") + "\t" + r"samtools idxstats " + i + " > " + sample_id + ".chr.txt\n")

    for i in os.listdir(cwd):
        if '.pmd3filter.bam' == i[-15:]:
            sample_id = i[:-19]
            print("\t> Reading distribution of reads per chromosome from '" + i + "'...")
            print("\t" + r"samtools idxstats " + i + " > " + sample_id + ".pmd.chr.txt")
            os.system(r"samtools idxstats " + i + " > " + sample_id + ".pmd.chr.txt")
            fout.write("\n\t> Reading distribution of reads per chromosome from '" + i + "'..." + time.strftime(
                "\n\t%H:%M:%S %Z ") + "\t" + r"samtools idxstats " + i + " > " + sample_id + ".pmd.chr.txt\n")

    print("\t> Distribution of reads per chromosome... Done!")
    fout.write("\n\t> Distribution of reads per chromosome... Done!\n\t-----\n")


def sexing():  ##Missing file_count error check
    """
    Performs sexing evaluation using Skoglund (2014) method. First the function creates and writes the authors' script into the working folder and then uses it as provided by them.

    Input:
    bar.q30.rmdup.bam

    Output:
    sexing_bar.txt

    """

    cwd = os.getcwd()
    script = open("XYkaryotyper.py", 'w')
    script.write(
        "#!/usr/bin/python2\n\n" + '"""' + "\nVersion:\t0.4\nAuthor:\t\tPontus Skoglund\nContact:\tpontus.skoglund@gmail.com\nDate:\t\t28 June 2013\n")
    script.write(
        "Citation:\tP. Skoglund, J.Stora, A. Gotherstrom, M. Jakobsson (2013)\n\t\tAccurate sex identification of ancient human remains using DNA shotgun sequencing.\n")
    script.write("\t\tJournal of Archaeological Science\n\n")
    script.write("Usage:\t\tpython ry_compute.py <SAM formatted data from stdin>\n\n")
    script.write("Example:\tsamtools view -q 30 mybamfile.bam | python XYkaryotyper.py\n\n")
    script.write(
        "\t\t(for specification on the SAM format and a the samtools suite, see Li, Handsaker et al. 2009, Bioinformatics)\n\n")
    script.write(
        "Output:\t\t[Total number of alignments in input] [Number of X and Y alignments identified] [R_y] [R_y standard error] [95% CI for R_Y] [Inferred sex]\n" + '"""' + "\n\n")
    script.write("import sys\nimport math\nfrom optparse import OptionParser\n\n")
    script.write('usage = "usage: %prog [options] <SAM formatted data from stdin>"\n')
    script.write("parser = OptionParser(usage=usage)\n")
    script.write(
        'parser.add_option("--chrXname", action="store", type="string", dest="chrXname",help="Identifier for the X chromosome in the SAM input (use if different than chrX, X etc)",default="X")\n')
    script.write(
        'parser.add_option("--chrYname", action="store", type="string", dest="chrYname",help="Identifier for the Y chromosome in the SAM input (use if different than chrY, Y etc)",default="Y")\n')
    script.write(
        'parser.add_option("--malelimit", action="store", type="float", dest="malelimit",help="Upper R_y limit for assignment as XY/male",default=0.075)\n')
    script.write(
        'parser.add_option("--femalelimit", action="store", type="float", dest="femalelimit",help="Lower R_y limit for assignment as XX/female",default=0.016)\n')
    script.write(
        'parser.add_option("--digits", action="store", type="int", dest="digits",help="Number of decimal digits in R_y output",default=4)\n')
    script.write(
        'parser.add_option("--noheader", action="store_true", dest="noheader",help="Do not print header describing the columns in the output",default=False)\n')
    script.write('(options, args) = parser.parse_args()\n\n')
    script.write("chrYcount=0\nchrXcount=0\ntotalcount=0\nfor line in sys.stdin:\n")
    script.write("\n\t#SAM header lines are skipped\n\tif line[0] == '@':continue\n\ttotalcount += 1\n\n")
    script.write("\tcol=line.split()\n\tchromosome=col[2]\n\n")
    script.write(
        "\tif options.chrYname in chromosome: chrYcount += 1\n\telif options.chrXname in chromosome: chrXcount += 1\n\n\n")
    script.write(
        "#compute R_y with 95% confidence interval\nn=chrYcount+chrXcount # total number of used alignments\nRy=1.0*chrYcount/n\nSE=math.sqrt((Ry*(1.0-Ry))/n)\nconfinterval=1.96*SE\n\n\n")
    script.write(
        "#use criteria to infer chromosomal sex\ngender='NA'\nif (Ry < options.femalelimit) and (Ry > options.malelimit):\n\tgender='Not Assigned'\nelif Ry==0.0:\n\tgender='consistent with XX'\n")
    script.write(
        "elif Ry+confinterval < options.femalelimit:\n\tgender='XX'\nelif Ry-confinterval > options.malelimit:\n\tgender='XY'\n")
    script.write(
        "elif Ry-confinterval > options.femalelimit and Ry+confinterval > options.malelimit:\n\tgender='consistent with XY but not XX'\nelif Ry-confinterval < options.femalelimit and Ry+confinterval < options.malelimit:\n")
    script.write(
        "\tgender='consistent with XX but not XY'\nelse:\n\tgender='Not Assigned'\n\nif options.noheader == False:\n\t" + r"print('Nseqs\tNchrY+NchrX\tNchrY\tR_y\tSE\t95% CI\tAssignment')" + "\n")
    script.write(
        r"print(str(totalcount)+'\t'+str(n)+'\t'+str(chrYcount)+'\t'+str(round(Ry,options.digits))+'\t'+str(round(SE,options.digits))+'\t'+str(round(Ry-confinterval,options.digits))+'-'+str(round(Ry+confinterval,options.digits))+'\t'+str(gender))")
    script.close()

    print("\n< SEXING - Skoglund's (2013) sex identification of ancient DNA human remains >")
    fout.write("\n< SEXING - Skoglund's (2013) sex identification of ancient DNA human remains >")

    for i in os.listdir(cwd):
        if '.rmdup.bam' == i[-10:]:
            sample_id, rest = re.split('\.q..',i)
            sep = re.findall('\.q..',i)[0]

            print("\t> Running sexing script for sample: " + sample_id)
            os.system('samtools view ' + i + ' | python XYkaryotyper.py > sexing_' + sample_id + '.txt')
            fout.write("\n\t> Running sexing script for sample: " + sample_id)
            fout.write(time.strftime(
                "\n\t%H:%M:%S %Z ") + '\tsamtools view ' + i + ' | python XYkaryotyper.py > sexing_' + sample_id + '.txt\n')

    print("\n\t> Sex identification of samples... Done!")
    fout.write("\n\t> Sex identification of samples... Done!\n\t-----\n")


def sexingEJ(single_out=True, st_out=True):  ##Missing file_count error check
    """
    Performs sexing evaluation using Jones (201X) method. First the function creates and writes the authors' script into the working folder and then uses it as provided by them.

    Input:
    bar.chr.txt

    Output:
    (Xy_all_samples.pdf)
    (Sexing_stats.txt)
    (bar.xy.pdf)

    """

    cwd = os.getcwd()
    print("\n< SEXING - Jones' (201X) sex identification of ancient DNA human remains >")
    fout.write("\n< SEXING - Jones' (201X) sex identification of ancient DNA human remains >")

    if single_out == True:
        arg1 = "single_out=TRUE"
    elif single_out == False:
        arg1 = "single_out=FALSE"
    else:
        print("\t> ERROR: Check function default arguments")
        quit()

    if st_out == True:
        arg2 = "st_out=TRUE"
    elif st_out == False:
        arg2 = "st_out=FALSE"
    else:
        print("\t> ERROR: Check function default arguments")
        quit()

    script = open("XYkaryotyper.R", 'w')
    script.write('setwd("')
    script.write(cwd)
    script.write('")')
    script.write(
        '\n\nargs=commandArgs(trailingOnly=TRUE)\n##Arguments order:\n##single_out=TRUE   st_out=TRUE\nargs=c("')
    script.write(arg1)
    script.write('","')
    script.write(arg2)
    script.write('")\n')
    script.write(
        '\n## Argument checking\nif(args[1]=="single_out=TRUE"){pdf("Xy_all_samples.pdf")}  ##For all graphs to be in same PDF\nif(args[2]=="st_out=TRUE"){final_tbl=data.frame(Sample=character(),XX_Prob=numeric(),XY_Prob=numeric(),p_value=numeric(),AssessXY=character(),Assessment=character(),stringsAsFactors=FALSE)}')
    script.write(
        '\n\nwdfiles=list.files(getwd())\nfor(i in wdfiles) {\n\tif(grepl(".chr.txt",i)==TRUE) {\n\t\t## Adapting input file to script and retrieving sample name\n\t\tname=strsplit(i,".chr.txt")')
    script.write(
        '\n\t\tname=name[[1]]\n\t\td=read.table(i,header=FALSE)\n\t\tnames(d)=c("chr","length",as.character(name))\n\t\td=`rownames<-`(d[-23, ], NULL)\n\t\td=d[-25,]\n\t\td=d[,-4]')
    script.write(
        '\n\t\tif(args[1]=="single_out=FALSE"){pdf(paste(name,".xy.pdf",sep=""))}  ##For individually separated graph PDF files\n\n\t\t###### START OF ORIGINAL SCRIPT ######')
    script.write(
        '\n\t\tcritt<-qt(0.05,df=20)\n\t\tprobs<-0\n\t\tratios<-0\n\n\t\tfor( i in 3:length(colnames(d)) )\n\t\t{\n\t\t\t# Start with plot')
    script.write(
        '\n\t\t\tplot(d[c(1:23),i]~d$length[c(1:23)],xlab="Chromsome length (Mb)",ylab="Number of reads", main=colnames(d)[i],pch=16,col="slategrey",xaxt="n")')
    script.write(
        '\n\t\t\tpoints(d$length[23],d[23,i],col="red",pch=16)\n\t\t\ttext(d$length[23],d[23,i],"X",col="red",pos=1)\n\t\t\ttext(d$length[1],d[1,i],"1",col="slategrey",pos=1, cex=0.5)')
    script.write(
        '\n\t\t\ttext(d$length[2],d[2,i],"2",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[3],d[3,i],"3",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[4],d[4,i],"4",col="slategrey",pos=1, cex=0.5)')
    script.write(
        '\n\t\t\ttext(d$length[5],d[5,i],"5",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[6],d[6,i],"6",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[7],d[7,i],"7",col="slategrey",pos=1, cex=0.5)')
    script.write(
        '\n\t\t\ttext(d$length[8],d[8,i],"8",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[9],d[9,i],"9",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[10],d[10,i],"10",col="slategrey",pos=1, cex=0.5)')
    script.write(
        '\n\t\t\ttext(d$length[11],d[11,i],"11",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[12],d[12,i],"12",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[13],d[13,i],"13",col="slategrey",pos=1, cex=0.5)')
    script.write(
        '\n\t\t\ttext(d$length[14],d[14,i],"14",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[15],d[15,i],"15",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[16],d[16,i],"16",col="slategrey",pos=1, cex=0.5)')
    script.write(
        '\n\t\t\ttext(d$length[17],d[17,i],"17",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[18],d[18,i],"18",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[19],d[19,i],"19",col="slategrey",pos=1, cex=0.5)')
    script.write(
        '\n\t\t\ttext(d$length[20],d[20,i],"20",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[21],d[21,i],"21",col="slategrey",pos=1, cex=0.5)\n\t\t\ttext(d$length[22],d[22,i],"22",col="slategrey",pos=1, cex=0.5)')
    script.write(
        '\n\t\t\tpoints(d$length[23],(2*d[23,i]),col="lightgrey",pch=16)\n\t\t\ttext(d$length[23],(2*d[23,i]),"2X",col="lightgrey",pos=3)\n\t\t\taxis(side=1,at=axTicks(1),labels=axTicks(1)/1000000)')
    script.write(
        '\n\n\t\t\t# Fit linear regression\n\t\t\tmodel<-lm(d[c(1:22),i]~d$length[c(1:22)])\n\t\t\tabline(model)\n\n\t\t\t# 95% CI lines for regression line')
    script.write(
        '\n\t\t\tprd <- predict( model, data.frame(d$length[c(1:22)]), interval="confidence", level=0.95)\n\t\t\tconfidence <- data.frame( cbind( d$length[c(1:22)], prd[,2], prd[,3] ) )')
    script.write(
        '\n\t\t\tconfordered <- confidence[order(confidence[,1]),]\n\t\t\tlines(confordered$X1,confordered$X2,lty=3)\n\t\t\tlines(confordered$X1,confordered$X3,lty=3)\n')
    script.write(
        '\n\t\t\t# Externally Studentized estimate of variance of residuals\n\t\t\tstudres<-rstudent(model)\n\t\t\tmodelwithX<-lm(d[c(1:23),i]~d$length[c(1:23)])\n\t\t\tstudreswithX<-rstudent(modelwithX)')
    script.write(
        '\n\t\t\tfemaleprob<-pt(studreswithX[23],df=20)\n\t\t\tmodelwithtwiceX<-lm(c(d[c(1:22),i],(d[23,i]*2))~d$length[c(1:23)])\n\t\t\tstudreswithtwiceX<-rstudent(modelwithtwiceX)')
    script.write(
        '\n\t\t\tmaleprob<-pt(studreswithtwiceX[23],df=20,lower.tail=F)\n\t\t\tratio<-femaleprob/maleprob\n\t\t\tprobs[i-2]<-femaleprob\n')
    script.write(
        '\n\t\t\t##Adition to script due to output table and blanks compatibility\n\t\t\tif (is.nan(ratio)==TRUE) { ratio[1]=0}\n\t\t\t##\n')
    script.write(
        '\n\t\t\tlm_eqn = function(df){\n\t\t\t\tmodel<-lm(d[c(1:22),i]~d$length[c(1:22)]);\n\t\t\t\teq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,list(a = format(coef(m)[1], digits = 2),b = format(coef(m)[2], digits = 2),r2 = format(summary(m)$r.squared, digits = 3)))')
    script.write(
        '\n\t\t\t\tas.character(as.expression(eq));\n\t\t\t}\n\n\t\t\t# Add descriptive legend\n\t\t\tif( ratio >= 1 )\n\t\t\t{\n\t\t\t\tlrt<--2*log(maleprob/femaleprob)        # <- stab at likelihood ratio test')
    script.write(
        '\n\t\t\t\tlrtp<-1-(pchisq(lrt,df=1))            # <- stab at computing probability based on LRT\n\t\t\t\t')
    script.write(
        'legend( "topleft",c(paste("Female (",signif(ratio,digits=3),"X more likely)",sep=""),paste("p = ", signif(lrtp,digits=3))), cex = 0.6)')
    script.write(
        '\n\t\t\t\tratios[i-2]<-ratio\n\t\t\t\t##Adition to original script due to output table\n\t\t\t\tif(args[2]=="st_out=TRUE"){\n\t\t\t\t\tassess=paste("Female (",signif(ratio,digits=3),"X more likely)",sep="")')
    script.write(
        '\n\t\t\t\t\tkar="XX"\n\t\t\t\t\tnewrow=data.frame(name,maleprob,femaleprob,signif(lrtp[[1]],digits=3),kar,assess,stringsAsFactors=FALSE)\n\t\t\t\t\tfinal_tbl=rbind(final_tbl,newrow)\n\t\t\t\t}\n\t\t\t\t##')
    script.write(
        '\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tlrt<--2*log(femaleprob/maleprob)        # <- stab at likelihood ratio test\n\t\t\t\tlrtp<-1-(pchisq(lrt,df=1))            # <- stab at computing probability based on LRT')
    script.write(
        '\n\t\t\t\tlegend( "topleft",c(paste("Male (",signif((1/ratio),digits=3),"X more likely)",sep=""),paste("p = ", signif(lrtp,digits=3))), cex = 0.6)')
    script.write(
        '\n\t\t\t\tratios[i-2]<-1/ratio\n\t\t\t\t##Adition to original script due to output table\n\t\t\t\tif(args[2]=="st_out=TRUE"){\n\t\t\t\t\tassess=paste("Male (",signif((1/ratio),digits=3),"X more likely)",sep="")')
    script.write(
        '\n\t\t\t\t\tkar="XY"\n\t\t\t\t\tnewrow=data.frame(name,maleprob,femaleprob,signif(lrtp[[1]],digits=3),kar,assess,stringsAsFactors=FALSE)\n\t\t\t\t\tfinal_tbl=rbind(final_tbl,newrow)\n\t\t\t\t}\n\t\t\t\t##\n\t\t\t}\n\t\t}')
    script.write(
        '\n\t\t#dev.off()\n\t\t######### END OF ORIGINAL SCRIPT #########\n\n\t\tif(args[1]=="single_out=FALSE"){dev.off()}  ##For individually separated graph PDF files\n\t}\n}\n\nif(args[2]=="st_out=TRUE"){\n\tempty_row=""')
    script.write(
        '\n\tcite_row=c("Based on script by Eppie Jones (201X).","","","","","")\n\tfinal_tbl=rbind(final_tbl,empty_row)\n\tfinal_tbl=rbind(final_tbl,cite_row)\n\tnames(final_tbl)=c("Sample","XX_Prob","XY_Prob","p_value","AssessXY","Assessment")\n\t')
    script.write(r'write.table(final_tbl,file="Sexing_stats.txt",sep="\t",quote=FALSE,row.names=FALSE)')
    script.write('\n}\n\ndev.off()')
    script.close()

    os.system("R CMD BATCH XYkaryotyper.R &")

    print("\n\t> Sex identification of samples... Done!")
    fout.write("\n\t> Sex identification of samples... Done!\n\t-----\n")


def mtDNAphymer(phyloTreelibraryPath, score1cutoff=0.04):
    """
    Uses the software PhyMer to estimate mitochondrial haplogroups and compiles the information from all samples into one tab-spaced file.

    Input:
    foo.q30.rmdup.bam

    Output:
    foo_mt.fastq
    foo.phymer.mtdna.txt
    mtDNA_Scores.txt

    :param phyloTreelibraryPath:
    (string)    Path to the PhyloTree to use

    :param score1cutoff:
    (string)    Lower threshold for accepting score1. Higher values give more precise results
    """

    ###################
    # ADD ARGUMENT TO CHECK PHYMER FOLDER INSTEAD OF PHYLOTREE, SO IT CAN BE USED FOR THE SNP LIST AS WELL
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< MT HAPLOGROUPS - Checking for mtDNA haplogroups using Phy-Mer >")
    fout.write("\n< MT HAPLOGROUPS - Checking for mtDNA haplogroups using Phy-Mer >")

    ## Checking if folder mtDNA, otherwise, create it
    if 'mtDNA' in os.listdir(cwd):
        pass
    else:
        newdir = cwd + '/' + "mtDNA"
        os.mkdir(newdir)
    cwdMT = cwd + "/mtDNA"

    ## New style of code - keep if no errors appear
    for i in os.listdir(cwd):
        if re.search("q30.rmdup.bam$", i) != None:
            sample_id, sep, rest = i.partition('.q30')
            print('\nsamtools view -b -h ' + i + ' chrM > ' + cwdMT + '/' + sample_id + '_mt.bam')
            os.system('samtools view -b -h ' + i + ' chrM > ' + cwdMT + '/' + sample_id + '_mt.bam')

    os.chdir(cwdMT)
    for i in os.listdir(cwdMT):
        if re.search("._mt.bam$", i) != None:
            # if '_mt.bam' == i[-7:]:        ## Using SAMTOOLS instead of picard for ease of use
            # sample_id, sep, rest = i.partition('_mt')
            # print('SamToFastq.jar INPUT='+i+' FASTQ='+sample_id+'_mt.fastq')
            # os.system('SamToFastq.jar INPUT='+i+' FASTQ='+sample_id+'_mt.fastq')
            sample_id, sep, rest = i.partition('_mt')
            print('samtools bam2fq ' + i + ' > ' + sample_id + '_mt.fastq')
            os.system('samtools bam2fq ' + i + ' > ' + sample_id + '_mt.fastq')

    for i in os.listdir(cwdMT):
        if re.search("_mt.fastq$", i) != None:
            # if '_mt.fastq' == i[-9:]:
            sample_id, sep, rest = i.partition('_mt')
            print(
                'Phy-Mer.py --verbose --print-ranking --def-snp=/home/dfernandes/Software/phy-mer-master/resources/Build_16_-_rCRS-based_haplogroup_motifs.csv ' + phyloTreelibraryPath + ' ' + i + ' > ' + sample_id + '.phymer.mtdna.txt')
            os.system(
                'Phy-Mer.py --verbose --print-ranking --def-snp=/home/dfernandes/Software/phy-mer-master/resources/Build_16_-_rCRS-based_haplogroup_motifs.csv ' + phyloTreelibraryPath + ' ' + i + ' > ' + sample_id + '.phymer.mtdna.txt')

    os.chdir(cwdMT)
    mtdnaout = open('mtDNA_Scores.txt', 'w')
    mtdnaout.write(
        "Sample\tHaplogroup1\tGeoLocation\tScore1\tHaplogroup2\tScore2\tHaplogroup3\tScore3\tHaplogroup4\tScore4\tHaplogroup5\tScore5")

    for i in os.listdir(cwdMT):
        if '.phymer.mtdna.txt' in i:
            sample_id, sep, rest = i.partition('.phymer')
            print("\n\t> Checking mtDNA haplogroup for sample " + sample_id)
            fout.write("\n\t> Checking mtDNA haplogroup for sample " + sample_id)
            phyfile = open(i, 'r')
            phyfiledata = phyfile.readlines()
            try:
                for line in phyfiledata:
                    if line[0] == "[":
                        first_index = phyfiledata.index(line)
                        break

                line1 = phyfiledata[first_index]
                line1 = line1.replace(' ', '')
                line1 = line1.replace("'", '')
                line1 = line1.replace("[", '')
                line1 = line1.replace("]", '')
                hap, sc1, sc2, sc3 = line1.split(',', 3)
                sc3 = sc3.rstrip("\n")
                if "\t" in sc3:
                    sc3, rest = sc3.split('\t', 1)
                line2 = phyfiledata[first_index + 1]
                line2 = line2.replace(' ', '')
                line2 = line2.replace("'", '')
                line2 = line2.replace("[", '')
                line2 = line2.replace("]", '')
                hap2, sc1b, sc2b, sc3b = line2.split(',', 3)
                sc3b = sc3b.rstrip("\n")
                if "\t" in sc3b:
                    sc3b, rest = sc3b.split('\t', 1)
                line3 = phyfiledata[first_index + 2]
                line3 = line3.replace(' ', '')
                line3 = line3.replace("'", '')
                line3 = line3.replace("[", '')
                line3 = line3.replace("]", '')
                hap3, sc1c, sc2c, sc3c = line3.split(',', 3)
                sc3c = sc3c.rstrip("\n")
                if "\t" in sc3c:
                    sc3c, rest = sc3c.split('\t', 1)
                line4 = phyfiledata[first_index + 3]
                line4 = line4.replace(' ', '')
                line4 = line4.replace("'", '')
                line4 = line4.replace("[", '')
                line4 = line4.replace("]", '')
                hap4, sc1d, sc2d, sc3d = line4.split(',', 3)
                sc3d = sc3d.rstrip("\n")
                if "\t" in sc3d:
                    sc3d, rest = sc3d.split('\t', 1)
                line5 = phyfiledata[first_index + 4]
                line5 = line5.replace(' ', '')
                line5 = line5.replace("'", '')
                line5 = line5.replace("[", '')
                line5 = line5.replace("]", '')
                hap5, sc1e, sc2e, sc3e = line5.split(',', 3)
                sc3e = sc3e.rstrip("\n")
                if "\t" in sc3e:
                    sc3e, rest = sc3e.split('\t', 1)
                if float(
                        sc1) < score1cutoff:  ## Checking if score1 is lower than defined cutoff. Too low values have no proper accuracy
                    pass
                else:
                    hapGeoList = open("/home/dfernandes/Software/PyScripts/mtDNAgeo.txt", 'r')
                    hapGeoListData = hapGeoList.readlines()
                    hapGeoDict = {}
                    for linha in hapGeoListData:
                        # print(linha)
                        hapGeo, sep, loc = linha.partition('\t ')
                        hapGeo = hapGeo.replace(" ", "")
                        loc = loc.rstrip('\t\n')
                        hapGeoDict[hapGeo] = loc
                    if hap in hapGeoDict:  ## Checking if main haplogroup has any geographic information available. If not, remove last letter/digit and check for the closest related haplogroup. Reduces up to a maximum of 4 letters/digits
                        # print("A: ",hap, hapGeoDict[hap])
                        print("\t\t" + hap)
                        fout.write("\n\t\t" + hap)
                        geoLoc = hapGeoDict[hap]
                    elif hap[:-1] in hapGeoDict:
                        print("\t\tFrom " + hap + " to " + hap[:-1])
                        fout.write("\n\t\tFrom " + hap + " to " + hap[:-1])
                        if hap[-1].isdigit() == True and hap[-2].isdigit() == True:
                            continue
                        else:
                            # print("B: ",hap[:-1], hapGeoDict[hap[:-1]])
                            geoLoc = hapGeoDict[hap[:-1]]
                    elif hap[:-2] in hapGeoDict:
                        print("\t\tFrom " + hap + " to " + hap[:-2])
                        fout.write("\n\t\tFrom " + hap + " to " + hap[:-2])
                        if hap[-2].isdigit() == True and hap[-3].isdigit() == True:
                            continue
                        else:
                            # print("C: ",hap[:-2], hapGeoDict[hap[:-2]])
                            geoLoc = hapGeoDict[hap[:-2]]
                    elif hap[:-3] in hapGeoDict:
                        print("\t\tFrom " + hap + " to " + hap[:-3])
                        fout.write("\n\t\tFrom " + hap + " to " + hap[:-3])
                        if hap[-3].isdigit() == True and hap[-4].isdigit() == True:
                            continue
                        else:
                            # print("D: ",hap[:-3], hapGeoDict[hap[:-3]])
                            geoLoc = hapGeoDict[hap[:-3]]
                    elif hap[:-4] in hapGeoDict:
                        print("\t\tFrom " + hap + " to " + hap[:-4])
                        fout.write("\n\t\tFrom " + hap + " to " + hap[:-4])
                        if hap[-4].isdigit() == True and hap[-5].isdigit() == True:
                            continue
                        else:
                            # print("E: ",hap[:-4], hapGeoDict[hap[:-4]])
                            geoLoc = hapGeoDict[hap[:-4]]

                    ## GEOLOC UNBOUNDLOCALERROR
                    try:
                        mtdnaout.write(
                            "\n" + sample_id + "\t" + hap + "\t" + geoLoc + "\t" + sc3 + "\t" + hap2 + "\t" + sc3b + "\t" + hap3 + "\t" + sc3c + "\t" + hap4 + "\t" + sc3d + "\t" + hap5 + "\t" + sc3e)
                    except UnboundLocalError:
                        pass

            except IndexError:
                pass


def Yfitter(index_file, pos_file, tree_file, mergeLanes = True):
    """
    Unfinished... But working!

    :param index_file:
    (string)    Path to reference indexed genome.

    :param tree_file:
    (string)    Path to mutation tree, provided with the software.

    :param mergeLanes:
    (boolean)   Merge info from reads from different lanes. NOTE: Current version does not work with "mergeLanes = False"!!

    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< Y HAPLOGROUPS - Checking for Y-chr haplogroups using Yfitter >")
    fout.write("\n< Y HAPLOGROUPS - Checking for Y-chr haplogroups using Yfitter >")

    if mapDamageRescaleSim == 1:
        fileExp = "mapdamaged.bam$"
    elif mapDamageRescaleSim == 0 and seqtkSim == 1:
        fileExp = "seqtk.bam$"
    elif mapDamageRescaleSim == 0 and seqtkSim == 0:
        fileExp = "q30.rmdup.bam$"

    ## Checking if any target file exists, otherwise, quit
    #check = any(fileExp == i[-14:] for i in os.listdir(cwd))
    check = any(re.search(fileExp,i[-len(fileExp):]) for i in os.listdir(cwd))
    if check == False:
        quit()
    if check == True:
        pass

    ## Checking if folder yfitter exists, otherwise, create it
    if 'yfitter' in os.listdir(cwd):
        pass
    else:
        newdir = cwd + '/' + "yfitter"
        os.mkdir(newdir)
    cwdY = cwd + "/yfitter"

    # NEW CODE!
    # DEBUGGING-Empty .bcf files might require the addition of flag "-A" to samtools mpielup command

    ## Extracts Y chromosome reads
    for i in os.listdir(cwd):
        if '.q30.rmdup.bam' == i[-14:]:
            sample_id, sep, rest = i.partition('.q30')
            print('\nsamtools view -b -h ' + i + ' chrY > ' + cwdY + '/' + sample_id + '_Y.bam')
            os.system('samtools view -b -h ' + i + ' chrY > ' + cwdY + '/' + sample_id + '_Y.bam')

    os.chdir(cwdY)
    for i in os.listdir(cwdY):
        if '_Y.bam' == i[-6:]:
            sample_id, sep, rest = i.partition('_Y.bam')
            ## Merge lanes if mergeLanes = True
            if mergeLanes == True:
                print('grep -v "^@RG" | sed "s/\\tRG:Z:[^\\t]*//" | samtools view - bo \n')
                print("samtools view -h " + cwdY + '/' + sample_id + '_Y.bam' + ' | grep -v "^@RG" | sed "s/\\tRG:Z:[^\\t]*//" | samtools view -bo ' + cwdY + '/' + sample_id + '_Ymer.bam')
                os.system("samtools view -h " + cwdY + '/' + sample_id + '_Y.bam' + ' | grep -v "^@RG" | sed "s/\\tRG:Z:[^\\t]*//" | samtools view -bo ' + cwdY + '/' + sample_id + '_Ymer.bam')
            if mergeLanes == True:
                samplePath = cwdY + '/' + sample_id + '_Ymer.bam'
            elif mergeLanes == False:
                samplePath = cwdY + '/' + sample_id + '_Y.bam'

            ## View/call all chrY genotypes
            print('\nsamtools mpileup -gf ' + index_file + ' ' + samplePath + ' > ' + cwdY + '/' + sample_id + '_Ys.bcf')
            os.system('samtools mpileup -gf ' + index_file + ' ' + samplePath + ' > ' + cwdY + '/' + sample_id + '_Ys.bcf')
            ## Index BCF file for call
            print('\nbcftools index ' + cwdY + '/' + sample_id + '_Ys.bcf')
            os.system('bcftools index ' + cwdY + '/' + sample_id + '_Ys.bcf')
            ## View/call defining genotypes
            ## NOTE: Currently using "view", but it might become useless in a few versions, being replaced by "call"
            print('\nbcftools view -R ' + pos_file + '.pos ' + cwdY + '/' + sample_id + '_Ys.bcf > ' + cwdY + '/' + sample_id + '_infosites.vcf')
            os.system('bcftools view -R ' + pos_file + '.pos ' + cwdY + '/' + sample_id + '_Ys.bcf > ' + cwdY + '/' + sample_id + '_infosites.vcf')
            ## Checks whether "infosites.vcf" file is of size 0, and if so skips next steps and use of Yfitter

            #REMOVE AFTER
            # sample_id = "M54A"
            # print(cwdY + '/' + sample_id + '_infosites.vcf')

            if os.path.getsize(cwdY + '/' + sample_id + '_infosites.vcf') != 0:
                ## Convert VCF to QCALL
                qcallOutPath = cwdY + '/' + sample_id + '_infosites.Qcall'
                qcallOut = open(qcallOutPath, "w")
                with open(cwdY + '/' + sample_id + '_infosites.vcf') as inVcf:
                    for line in inVcf:
                        if line[0] == "#":
                            pass
                        elif line[0] != "#":
                            #print(line)
                            #print(line.split()[3])
                            #print(line.split()[4])
                            #print(line.split()[5])
                            if "INDEL" in line:
                                pass
                            elif "DP=" in line:
                                bases = line.split()[3:5]
                                #print(bases)
                                #print(len(bases))
                                #print("AAA", bases[len(bases) - 1].split(","))
                                #print("AAA", bases[len(bases) - 1].split(",")[0])
                                if "," in bases[1]:
                                    #print("BBB", bases[len(bases) - 1].split(",")[0])
                                    bases[len(bases) - 1] = bases[len(bases) - 1].split(",")[0]
                                if bases[len(bases)-1] == "*" or bases[len(bases)-1] == "<*>":
                                    bases = bases[0:len(bases)-1]
                                if len(bases[0]) > 1:
                                    pass
                                elif len(bases[0]) == 1:
                                    #print(bases)
                                    outLine = []
                                    outLine.append(line.split()[0])
                                    outLine.append(line.split()[1])
                                    outLine.append(line.split()[3])
                                    outLine.append(line.split()[7].split(";")[0].split("DP=")[1])
                                    outLine.append("37")
                                    outLine.append("0")
                                    #print(outLine)
                                    likes = line.split()[9].split(",")
                                    #print(likes)
                                    if len(bases) == 1:
                                        #print("PartA")
                                        if bases[0] == "A":
                                            outLine.append(likes[0])
                                            outLine.append(likes[1])
                                            outLine.append(likes[1])
                                            outLine.append(likes[1])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                        if bases[0] == "C":
                                            outLine.append(likes[2])
                                            outLine.append(likes[1])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[0])
                                            outLine.append(likes[1])
                                            outLine.append(likes[1])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                        if bases[0] == "G":
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[1])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[1])
                                            outLine.append(likes[2])
                                            outLine.append(likes[0])
                                            outLine.append(likes[1])
                                            outLine.append(likes[2])
                                        if bases[0] == "T":
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[1])
                                            outLine.append(likes[2])
                                            outLine.append(likes[2])
                                            outLine.append(likes[1])
                                            outLine.append(likes[2])
                                            outLine.append(likes[1])
                                            outLine.append(likes[0])
                                    elif len(bases) == 2:
                                        #print("PartB")
                                        #likes = likes[0:3]
                                        #print(likes)
                                        if bases[0] == "A":
                                            if bases[1] == "C":
                                                outLine.append(likes[0])
                                                outLine.append(likes[1])
                                                outLine.append(likes[3])
                                                outLine.append(likes[3])
                                                outLine.append(likes[2])
                                                outLine.append(likes[4])
                                                outLine.append(likes[4])
                                                outLine.append(likes[5])
                                                outLine.append(likes[5])
                                                outLine.append(likes[5])
                                            if bases[1] == "G":
                                                outLine.append(likes[0])
                                                outLine.append(likes[3])
                                                outLine.append(likes[1])
                                                outLine.append(likes[3])
                                                outLine.append(likes[5])
                                                outLine.append(likes[4])
                                                outLine.append(likes[5])
                                                outLine.append(likes[2])
                                                outLine.append(likes[4])
                                                outLine.append(likes[5])
                                            if bases[1] == "T":
                                                outLine.append(likes[0])
                                                outLine.append(likes[3])
                                                outLine.append(likes[3])
                                                outLine.append(likes[1])
                                                outLine.append(likes[5])
                                                outLine.append(likes[5])
                                                outLine.append(likes[4])
                                                outLine.append(likes[5])
                                                outLine.append(likes[4])
                                                outLine.append(likes[2])
                                        if bases[0] == "C":
                                            if bases[1] == "A":
                                                outLine.append(likes[2])
                                                outLine.append(likes[1])
                                                outLine.append(likes[4])
                                                outLine.append(likes[4])
                                                outLine.append(likes[0])
                                                outLine.append(likes[3])
                                                outLine.append(likes[3])
                                                outLine.append(likes[5])
                                                outLine.append(likes[5])
                                                outLine.append(likes[5])
                                            if bases[1] == "G":
                                                outLine.append(likes[5])
                                                outLine.append(likes[3])
                                                outLine.append(likes[4])
                                                outLine.append(likes[5])
                                                outLine.append(likes[0])
                                                outLine.append(likes[1])
                                                outLine.append(likes[3])
                                                outLine.append(likes[2])
                                                outLine.append(likes[4])
                                                outLine.append(likes[5])
                                            if bases[1] == "T":
                                                outLine.append(likes[5])
                                                outLine.append(likes[3])
                                                outLine.append(likes[5])
                                                outLine.append(likes[4])
                                                outLine.append(likes[0])
                                                outLine.append(likes[3])
                                                outLine.append(likes[1])
                                                outLine.append(likes[5])
                                                outLine.append(likes[4])
                                                outLine.append(likes[2])
                                        if bases[0] == "G":
                                            if bases[1] == "A":
                                                outLine.append(likes[2])
                                                outLine.append(likes[4])
                                                outLine.append(likes[1])
                                                outLine.append(likes[4])
                                                outLine.append(likes[5])
                                                outLine.append(likes[3])
                                                outLine.append(likes[5])
                                                outLine.append(likes[0])
                                                outLine.append(likes[3])
                                                outLine.append(likes[5])
                                            if bases[1] == "C":
                                                outLine.append(likes[5])
                                                outLine.append(likes[4])
                                                outLine.append(likes[3])
                                                outLine.append(likes[5])
                                                outLine.append(likes[2])
                                                outLine.append(likes[1])
                                                outLine.append(likes[4])
                                                outLine.append(likes[0])
                                                outLine.append(likes[3])
                                                outLine.append(likes[5])
                                            if bases[1] == "T":
                                                outLine.append(likes[5])
                                                outLine.append(likes[5])
                                                outLine.append(likes[3])
                                                outLine.append(likes[4])
                                                outLine.append(likes[5])
                                                outLine.append(likes[3])
                                                outLine.append(likes[4])
                                                outLine.append(likes[0])
                                                outLine.append(likes[1])
                                                outLine.append(likes[2])
                                        if bases[0] == "T":
                                            if bases[1] == "A":
                                                outLine.append(likes[2])
                                                outLine.append(likes[4])
                                                outLine.append(likes[4])
                                                outLine.append(likes[1])
                                                outLine.append(likes[5])
                                                outLine.append(likes[5])
                                                outLine.append(likes[3])
                                                outLine.append(likes[5])
                                                outLine.append(likes[3])
                                                outLine.append(likes[0])
                                            if bases[1] == "C":
                                                outLine.append(likes[5])
                                                outLine.append(likes[4])
                                                outLine.append(likes[5])
                                                outLine.append(likes[3])
                                                outLine.append(likes[2])
                                                outLine.append(likes[4])
                                                outLine.append(likes[1])
                                                outLine.append(likes[5])
                                                outLine.append(likes[3])
                                                outLine.append(likes[0])
                                            if bases[1] == "G":
                                                outLine.append(likes[5])
                                                outLine.append(likes[5])
                                                outLine.append(likes[4])
                                                outLine.append(likes[3])
                                                outLine.append(likes[5])
                                                outLine.append(likes[4])
                                                outLine.append(likes[3])
                                                outLine.append(likes[2])
                                                outLine.append(likes[1])
                                                outLine.append(likes[0])
                                    outLine.append(sample_id)
                                    qcallOut.write("\t".join(outLine))
                                    qcallOut.write("\n")
                qcallOut.close()
                ## Run Yfitter
                print('\nYfitter ' + tree_file + '.xml ' + cwdY + '/' + sample_id + '_infosites.Qcall > ' + cwdY + '/' + sample_id + '_Yfitter')
                os.system('Yfitter ' + tree_file + '.xml ' + cwdY + '/' + sample_id + '_infosites.Qcall > ' + cwdY + '/' + sample_id + '_Yfitter')


def yleaf(refPath, refVersion, minReadThres = 1, minQualThres = 20, minPercBaseAgreementThres = 90, numThreads = 1):
    """
    Predicts Y chromosome haplogroups using Yleaf software (https://github.com/dmontielg/Yleaf).

    :input:
    foo.qxx.rmdup.bam

    :param refPath:
    (string)    Path to reference indexed genome.

    :param refVersion:
    (string)    Version of the reference genome to use (accepts "hg19", "hg38").

    :param minReadThres:
    (int)   Minimum number of reads per position (default = 1).

    :param minQualThres:
    (int)   Minimum accepted base quality (accepts values between 10 and 39) (default = 20).

    :param minPercBaseAgreementThres:
    (int)   Minimum percentage of allele agreement per position (default = 90).

    :param numThreads:
    (int)   Number of CPU threads to use (default = 1).

    :return:
    Yleaf (folder)
        foo (folder)
            foo.hg
            foo.qxx.rmdup (folder)
                foo.qxx.rmdup.chr
                foo.qxx.rmdup.fmf
                foo.qxx.rmdup.log
                foo.qxx.rmdup.out
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< Y LEAF - Estimating Y chromosome haplogroups >")
    fout.write("\n< Y LEAF - Estimating Y chromosome haplogroups >")

    newdir = cwd + '/' + "yleaf"
    if 'yleaf' in os.listdir(cwd):
        shutil.rmtree(newdir)

    os.mkdir(newdir)
    cwdY = cwd + "/yleaf"

    ## Find Yleaf's working folder and define position file (working for Linux Mint)
    positionFileLoc = os.popen("whereis Yleaf").read().split(" ")[1].rsplit(r"/",1)[0] + "/Position_files/"+refVersion+".txt"

    ## Main loop
    for i in os.listdir(cwd):
        if re.search('\.q\d\d\.rmdup\.bam$', i) != None:
            print(i)
            sample_id, sep, rest = i.partition('.q30')
            ## Always delete existing Yleaf sample folder in order to avoid software prompt
            if sample_id +'_yleaf' in os.listdir(cwd):
                shutil.rmtree(sample_id +'_yleaf')
                input("Removed folder")
            print("Yleaf.py -bam " + i + " -f " + refPath + " -pos " + positionFileLoc + " -out " + sample_id + "_yleaf -r " + str(minReadThres) + " -q " + str(minQualThres) + " -b " + str(minPercBaseAgreementThres) + " -t " + str(numThreads))
            os.system("Yleaf.py -bam " + i + " -f " + refPath + " -pos " + positionFileLoc + " -out " + sample_id + "_yleaf -r " + str(minReadThres) + " -q " + str(minQualThres) + " -b " + str(minPercBaseAgreementThres) + " -t " + str(numThreads))
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + '\tYleaf.py -bam ' + i + " -f " + refPath + " -pos " + positionFileLoc + " -out " + sample_id + "_yleaf -r " + str(minReadThres) + " -q " + str(minQualThres) + " -b " + str(minPercBaseAgreementThres) + " -t " + str(numThreads) + '\n')

    ## Move folders into "yleaf" folder
    for i in os.listdir(cwd):
        if re.search('_yleaf$', i) != None:
            os.system('mv ' + i + ' ' + cwdY)

    ## Exclude C>T and G>A sites and generate report (in R)
    script = open("Yleaf_Batch_Read_Exclude_CTAG.R", 'w')
    script.write('setwd("')
    script.write(cwdY)
    script.write('")\n\n')
    script.write("summary_df = data.frame(Sample=character(),Ychr_Valid_Markers=character(),Ychr_Valid_Markers_NoCT_NoAG=character(),Ychr_Valid_Derived_Markers_NoCT_NoAG=character(),Most_Downstream_Haplogroup=character(),Total_Derived_Markers_For_Upstream_Major_Haplogroup=character(),stringsAsFactors=FALSE)\n")
    script.write("\nfor(i in list.dirs(recursive = F)) {\n")
    script.write("\tsubFolder = list.dirs(paste0(getwd(),substring(i,2)))[2]\n\t")
    script.write(r'sampleID = strsplit(subFolder,"/")[[1]][length(strsplit(subFolder,"/")[[1]])]')
    script.write("\n\t")
    script.write(r'if(file.exists(paste0(subFolder,"/",sampleID,".out"))) {')
    script.write("\n\t\t")
    script.write(r'readIn = read.table(file = paste0(subFolder,"/",sampleID,".out"),header = T,stringsAsFactors = F)')
    script.write('\n\t\tvalidMarkers = length(readIn$marker_name)\n\t\tcts = which(readIn$mutation == "C->T")\n\t\tgas = which(readIn$mutation == "G->A")\n\t\treadIn = readIn[-c(cts,gas),]\n\t\t')
    script.write(r'validMarkersNoCTAG = length(readIn$marker_name)')
    script.write('\n\t\t')
    script.write(r'deriveds = readIn[which(readIn$state == "D"),]')
    script.write('\n\t\t')
    script.write(r'validDerivedMarkersNoCTAG = length(deriveds$marker_name)')
    script.write('\n\t\t')
    script.write(r'mostDownstreamHap = deriveds$haplogroup[length(deriveds$haplogroup)]')
    script.write('\n\t\t')
    script.write(r'otherSameMajorDeriveds = length(which(substr(mostDownstreamHap,start = 1,stop = 1) == substr(deriveds$haplogroup,start = 1,stop=1)))')
    script.write('\n\t\tsample_row=c(sampleID,validMarkers,validMarkersNoCTAG,validDerivedMarkersNoCTAG,mostDownstreamHap,paste0(otherSameMajorDeriveds," (",substr(mostDownstreamHap,start = 1,stop = 1),")"))\n\t\t')
    script.write('summary_df=rbind(summary_df,sample_row,stringsAsFactors = FALSE)\n\t\t')
    script.write(r'write.table(x = deriveds, file = paste0(subFolder,"/",sampleID,".out_derived_NoCTAG"), quote = F, sep = "\t",row.names = F,col.names = T)')
    script.write('\n\t} else {\n\t\t')
    script.write(r'sample_row=c(sampleID,NA,NA,NA,NA,paste0(NA," (",NA,")"))')
    script.write('\n\t\tsummary_df=rbind(summary_df,sample_row,stringsAsFactors = FALSE)\n\t}\n}\ncolnames(summary_df) = c("Sample","Ychr_Valid_Markers","Ychr_Valid_Markers_NoCT_NoAG","Ychr_Valid_Derived_Markers_NoCT_NoAG","Most_Downstream_Haplogroup","Total_Derived_Markers_(For_Upstream_Major_Haplogroup)")\n')
    script.write(r'write.table(x = summary_df,file = paste0(paste0(strsplit(getwd(),"/")[[1]][1:length(strsplit(getwd(),"/")[[1]])-1],collapse = "/"), "/Summary_Yleaf.txt"), quote = F, col.names = T, row.names = F, sep="\t")')
    script.close()

    os.system("R CMD BATCH Yleaf_Batch_Read_Exclude_CTAG.R")
    fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tR CMD BATCH Yleaf_Batch_Read_Exclude_CTAG.R")
    os.remove(cwd + "/Yleaf_Batch_Read_Exclude_CTAG.R")
    os.remove(cwd + "/Yleaf_Batch_Read_Exclude_CTAG.Rout")








def readlengths():  ##Missing file_count error check
    """
    (path) -> various_files

    Analyzes average read length and standard deviation on '*q30.rmdup.bam' files using PICARDTOOLS and compiles the batched info with an R script.
    Results output to "Reads_Lengths.txt".
    Requires 'jarwrapper'.


    Input:
    foo.q30.rmdup.bam

    Output:
    foo.readlens.txt
    Reads_Lengths.txt

    Example:
    >>> readlengths()
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< READ_LENGTHS - Average read lengths >")
    fout.write("\n< READ_LENGTHS - Average read lengths >")

    for i in os.listdir(cwd):
        if '.rmdup.bam' == i[-10:]:
            #sample_id, sep, rest = i.partition('.q30')
            sample_id, rest = re.split('\.q..', i)
            print("\t> Running read length analysis for sample: " + sample_id)
            os.system('picard.jar CollectQualityYieldMetrics ' + ' I=' + i + ' O=' + sample_id + '.readlens.txt')
            print('picard.jar CollectQualityYieldMetrics ' + ' I= ' + i + ' O=' + sample_id + '.readlens.txt')
            fout.write("\t> Running read length analysis for sample: " + sample_id)
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + '\tpicard.jar CollectQualityYieldMetrics ' + ' I= ' + i + ' O=' + sample_id + '.readlens.txt\n')

    script = open("read_lengths_batch.R", 'w')
    script.write('getwd()\nsetwd("')
    script.write(cwd)
    script.write('")')
    script.write('\nwdfiles=list.files(getwd())')
    script.write("\n\nfinal_df = data.frame(Sample=character(),Mean_Length=character(),Standard_Deviation=character(),stringsAsFactors=FALSE)")
    script.write('\n\nfor(i in wdfiles) {\n\tif(grepl("readlens.txt$",i)==TRUE) {\n\t\tfil=file(i,open="rt")\n\t\tres=readLines(fil)')
    script.write('\n\n\t\t## Get sample name\n\t\tsample_id = strsplit(strsplit(strsplit(res[2]," ")[[1]][3],"=")[[1]][2],".q...rmdup.bam")[[1]][1]')
    script.write('\n\n\t\t## Get lengths\n\t\tmean_len = strsplit(res[grep(pattern = "^TOTAL_READS", res)+1],')
    script.write(r'"\\\t")[[1]][3]')
    script.write('\n\t\tsample_row=c(sample_id,mean_len)\n\t\tfinal_df=rbind(final_df,sample_row,stringsAsFactors = FALSE)\n\t\tclose(fil)\n\t}\n}\n\ncolnames(final_df) = c("Sample","Mean_Length")\n\n')
    script.write(r'write.table(final_df,file="Read_Lenghts.txt",sep="\t",quote=FALSE,row.names=FALSE)')
    script.close()

    os.system("R CMD BATCH read_lengths_batch.R &")

    print("\n\t> Read length analysis... Done!")
    fout.write("\n\t> Read length analysis... Done!\n\t-----\n")


def pca_lc(filePath, datasetFolder, mapFileDatasetPath, intervalListFile, indexGen=False, makepileup=True, minNumReadsPileup=15000, makeped=True,
           pedfamily=False, makemap=True, conv2bin=True, mergeIt=True, reduceDSsnp=False, genoPlink=1,
           reduceDSind=False, popList=None, PCAparamFile=True, plot=True):
    """
    This pipeline processes input *.bam files to do PCA analysis using the 'smartpca' software. First, uses GATK to call the SNPs present in the given dataset.
    Then it converts this into plink data, which is used to merge with the refernce dataset and proceed with the analysis.
    Please carefully read the description of each parameter if you have any doubts or errors.

    Input:
    Foo.q30.rmdup.bam
    intervalListFile

    Output:
    Foo.pileupgatk.haak2015.txt
    Foo.ped
    Foo.map
    Foo.bed/bim/fam
    Foomergedex.bed/bim/fam
    Foomergedex.eigenstratgeno/ind/snp
    Foomergedex.evec/eval

    :param filePath:
    (str)   File path to index genome used for the samples and datasets. To be used for new indexing step needed by GATK.

    :param datasetFolder:
    (str)   Path to the folder containing the dataset files. Note this is not the dataset name itself, just the folder.

    :param mapFileDatasetPath:
    (str)   Path to the *.bim file from the selected dataset.

    :param intervalListFile:
    (str)   Path to interval file listing the SNP positions of the dataset, to be used for creation of pileup file with GATK.

    :param indexGen:
    (boolean)   Index source genome with samtools faidx?

    :param makepileup:
    (boolean)   Create pileup files with GATK?

    :param minNumReadsPileup:
    (int)   Lower threshold of read numbers to allow Pileup

    :param makeped:
    (boolean)   Make PED files?

    :param pedfamily:
    (boolean)
    (str)   Include PED family information for the samples analysed? Be aware all samples will have the same family ID. Must be a string, or set to pedFamily=False if not used.

    :param makemap:
    (boolean)   Make MAP files?

    :param conv2bin:
    (boolean)   Convert PED+MAP to binary BED+BIM+FAM?

    :param mergeIt:
    (boolean)   Merge dataset with sample?

    :param reduceDSsnp:
    (boolean)   Reduce dataset's SNPs to match those of your sample? Uses R to create a list of SNPs from the sample, and applies it to plink, reducing the dataset.

    :param genoPlink:
    (int)   Apply missing genotype rate SNP exclusion step. Default is 1, which includes all SNPs.

    :param reduceDSind:
    (boolean)
    (str)   Use dataset with reduced populations/individuals? Must be a string with the name of the dataset (without file extenstions) if it is to be used, or set to reduceDSind=False if not used.
    CAUTION: "False" parameter might broken code for now, always provide dataset to reduce

    :param popList:
    (str)   Path to population list for creation of PCA using smartpca. Default is popList=None.

    :param PCAparamFile:
    (boolean)   Create PCA parameter file?

    :param plot:
    (boolean)   Plot the PCAs?

    :return:
    """

    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< ROADTOPCA - Pipeline for creation of PCA plots from bam files >")
    fout.write("\n< ROADTOPCA - Pipeline for creation of PCA plots from bam files >")
    cwdgen, sep, id = filePath.rpartition('/')

    if indexGen == True:
        cwdgen, sep, id = filePath.rpartition('/')
        print("\t> Indexing reference genome " + "'" + id + "'" + " in '" + cwd + "'...")
        os.system('samtools faidx ' + filePath)
    elif indexGen == False:
        pass

    dataset = datasetFolder.rsplit("/", 1)
    dataset = dataset[1]
    print("\t> Dataset used:", dataset)

    mapFile = mapFileDatasetPath.rsplit("/", 1)
    mapFile = mapFile[1]


    ## PILEUP conversion w/ GATK
    if makepileup == True:
        for i in os.listdir(cwd):
            if '.rmdup.bam' == i[-10:]:
                #sample_id, sep, rest = i.partition('.q30')
                sample_id, rest = re.split('\.q..', i)
                bamlength = int(os.popen('samtools view -c ' + i).read())
                print("\n\t> Collecting PILEUP data from sample: " + sample_id)
                if i + '.bai' in os.listdir(cwd):
                    if bamlength > minNumReadsPileup:
                        print(
                            '\tgatk \n\tPileup \n\t-R ' + filePath + ' \n\t-L ' + intervalListFile + ' \n\t-I ' + i + ' \n\t-O ' + sample_id + '.pileupgatk.' + dataset + '.txt')
                        os.system(
                            'gatk Pileup -R ' + filePath + ' -L ' + intervalListFile + ' -I ' + i + ' -O ' + sample_id + '.pileupgatk.' + dataset + '.txt')
                        fout.write("\n\t> Collecting PILEUP data from sample: " + sample_id + time.strftime(
                            "\n\t%H:%M:%S %Z ") + '\tgatk \n\t\t\tPileup \n\t\t\t-R ' + filePath + ' \n\t\t\t-L ' + intervalListFile + ' \n\t\t\t-I ' + i + ' \n\t\t\t-O ' + sample_id + '.pileupgatk.' + dataset + '.txt\n')
                    else:
                        fout.write("\n\t> Collecting PILEUP data from sample: " + sample_id)
                        print(
                            "\t\t> WARNING: Input file has less than 15000 reads. \n\t\t> Expected number of SNPs might not be informative. Change 'minNumReadsPileup' parameter if you want to continue.\n\t\t> Skipping PCA analysis for sample " + sample_id)
                        fout.write(
                            "\n\t\t> WARNING: Input file has less than 15000 reads. \n\t\t> Expected number of SNPs might not be informative. Change 'minNumReadsPileup' parameter if you want to continue.\n\t\t> Skipping PCA analysis for sample " + sample_id + "\n")
                        pass
                else:
                    print(
                        "\t> ERROR: Index file (*.bai) for bam file '" + i + "' not found. Aborting pileup creation...")
                    fout.write(
                        "\t> ERROR: Index file (*.bai) for bam file '" + i + "' not found. Aborting pileup creation...")
                    pass
    elif makepileup == False:
        pass
    elif makepileup != True and makepileup != False:
        print("\t\t> ERROR: 'makepileup' argument has invalid input. Must be 'True' or 'False'.")

    ## PED creation from PILEUP GATK
    if makeped == True:
        dataset = datasetFolder.rsplit("/", 1)
        dataset = dataset[1]
        for i in os.listdir(cwd):
            if '.rmdup.bam' == i[-10:]:
                #sample_id, sep, rest = i.partition('.q30')
                sample_id, rest = re.split('\.q..', i)
                found_file = 0
                print("\t> Creating PED file from sample: " + sample_id)
                fout.write("\n\t> Creating PED file from sample: " + sample_id)
                for i in os.listdir(cwd):
                    file_wanted = sample_id + '.pileupgatk.' + dataset + ".txt"
                    if file_wanted == i:
                        # if sample_id+'.pileupgatk.' in i and dataset in i:
                        found_file = 1
                        sample_id, sep, rest = i.partition('.pileup')
                        pedout = open(sample_id + '.ped', 'w')
                        pileupin = open(i, 'r')
                        snpcount = sum(1 for _ in pileupin)
                        print("\tSource pileup file:", i)
                        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tSource pileup file:" + i)
                        print("\tDetected SNPs:", snpcount)
                        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tDetected SNPs:" + str(snpcount))
                        if pedfamily == False:  ## Writing all samples as MALES by default, as no-sex allows for some complications to occur downstream
                            #                            print("\tFamily ID, Sample ID, Paternal and Maternal IDs, Sex, Affection: "+sample_id+' '+sample_id+' 0 0 1 1 ') ## 1 is MALE, 2 is FEMALE, 0 is UNKNOWN
                            print(
                                "\tFamily ID, Sample ID, Paternal and Maternal IDs, Sex, Affection: " + sample_id + ' ' + sample_id + ' 0 0 1 1 ')
                            fout.write(time.strftime(
                                "\n\t%H:%M:%S %Z ") + "\tFamily ID, Sample ID, Paternal and Maternal IDs, Sex, Affection: " + sample_id + ' ' + sample_id + ' 0 0 1 1 \n')
                            pedout.write(sample_id + ' ' + sample_id + ' 0 0 1 1 ')
                        else:
                            print(
                                "\tFamily ID, Sample ID, Paternal and Maternal IDs, Sex, Affection: " + pedfamily + ' ' + sample_id + ' 0 0 1 1 ')
                            pedout.write(pedfamily + ' ' + sample_id + ' 0 0 1 1 ')
                        pileupin = open(i, 'r')
                        counter = 1
                        countbads = 0
                        for line in pileupin:
                            if counter <= snpcount:
                                line = line.rstrip()
                                chro, coord, refb, rbase, qual = line.split(' ')
                                counter += 1
                                ## Identify if there is more than one allele for this position
                                if len(rbase) > 1:
                                    listBases = list(rbase)
                                    listQual = list(qual)
                                    listPhred = []
                                    ## Only get alleles with phred quality > 30
                                    for i in listQual:
                                        phred = ord(i) - 33
                                        listPhred.append(phred)
                                    listOfGoodBases = []
                                    counterPhred = 1
                                    while counterPhred <= len(listPhred):
                                        for ph in listPhred:
                                            if ph >= 30:
                                                listOfGoodBases.append(listBases[counterPhred - 1])
                                                counterPhred += 1
                                            else:
                                                counterPhred += 1
                                                pass
                                    if len(listOfGoodBases) != 0:
                                        chosenbase = random.choice(listOfGoodBases)
                                        pedout.write(chosenbase + ' ' + chosenbase + ' ')
                                    else:
                                        # input("B - This position did not have any good quality base")
                                        countbads += 1
                                        pedout.write('0 0 ')
                                else:
                                    rbasePhred = ord(qual) - 33
                                    if rbasePhred >= 30:
                                        pedout.write(rbase + ' ' + rbase + ' ')
                                    else:
                                        # input("C - This position has only 1 read and it is not good")
                                        pedout.write('0 0 ')
                                        countbads += 1
                            else:  ## When reaching the last line, pass
                                pass
                        pedout.close()
                        print("\tSNPs skipped due to low quality: " + str(countbads) + "\n")
                        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tSNPs skipped due to low quality: " + str(
                            countbads) + "\n")
                if found_file != 1:
                    print(
                        "\t\t> ERROR: Pileup file for sample '" + sample_id + "' not found. Aborting PED file creation...\n")
                    fout.write(
                        "\n\t\t> ERROR: Pileup file for sample '" + sample_id + "' not found. Aborting PED file creation...\n")
    elif makeped == False:
        pass
    elif makeped != True and makeped != False:
        print("\t> ERROR: 'makeped' argument has invalid input. Must be 'True' or 'False'.")


## Create MAP file
    if makemap == True:
        listOfSamples = []
        for i in os.listdir(cwd):
            if ".ped" in i:
                sample_id, sep, rest = i.partition('.ped')
                for i in os.listdir(cwd):
                    file_wanted_pile = sample_id + '.pileupgatk.' + dataset + ".txt"
                    if file_wanted_pile == i:
                        listOfSamples.append(sample_id)

        ##Creating MAP files with R (20X faster than dictionary method)
        print("\t> Creating MAP files")
        fout.write("\n\n\t> Creating MAP files")
        script = open("create_maps.R", 'w')
        script.write('setwd("')
        script.write(cwd)
        script.write('")\nsampleVec = c(')
        counter = 1
        for samp in listOfSamples:
            if counter == 1:
                script.write('"')
                script.write(samp)
                script.write('"')
            else:
                script.write(',"')
                script.write(samp)
                script.write('"')
            counter += 1
        script.write(')\n\nloadOriginalMap = read.table("')
        script.write(mapFileDatasetPath)
        script.write('", stringsAsFactors = F)')
        script.write('\n\nfor (i in sampleVec) {\n\tloadPileUp = read.table(paste0(i, ".pileupgatk.')
        script.write(dataset)
        script.write('.txt"), stringsAsFactors = F, fill = T, quote = "")\n\t')
        script.write('\n\n\tloadPileUp$V6 = paste0(gsub("[a-z]", "", loadPileUp$V1), "-", loadPileUp$V2)\n\tif(length(grep("X",loadPileUp$V6)) >0){\n\t\tloadPileUp$V6 = gsub("X","23",loadPileUp$V6)\n\t}') 
        script.write('\n\tif(length(grep("Y",loadPileUp$V6)) >0){\n\t\tloadPileUp$V6 = gsub("Y","24",loadPileUp$V6)\n\t}\n\tloadOriginalMap$V7= paste0(loadOriginalMap$V1, "-", loadOriginalMap$V4)')
        script.write('\n\n\tsaveMap = loadOriginalMap[which(loadOriginalMap$V7 %in% loadPileUp$V6), ]\n\tsaveMap$V5 = NULL\n\tsaveMap$V6 = NULL\n\tsaveMap$V7 = NULL')
        script.write('\n\n\twrite.table(saveMap, paste0(i, ".map"), quote = F, sep = "\t", row.names = F, col.names = F)\n}')
        script.close()
        os.system("R CMD BATCH create_maps.R")
        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tR CMD BATCH create_maps.R")
        os.remove(cwd+"/create_maps.R")
        os.remove(cwd + "/create_maps.Rout")

    elif makemap == False:
        pass
    elif makemap != True and makemap != False:
        print("\t> ERROR: 'makemap' argument has invalid input. Must be 'True' or 'False'.")
        fout.write("\n\t> ERROR: 'makemap' argument has invalid input. Must be 'True' or 'False'.")

        ## MAP IS FINE, FINISH WITH FOUT.WRITE AND TIME STAMPS AND MAYBE PROGRESS BAR

    ## Convert to binary
    if conv2bin == True:
        if "bin4pca" in os.listdir(cwd):
            pass
        else:
            newdir = cwd + '/' + "bin4pca"
            os.mkdir(newdir)
        for i in os.listdir(cwd):
            if ".ped" == i[-4:]:
                sample_id, sep, rest = i.partition('.ped')
                print("\t> Converting plink data to binary format for " + sample_id)
                fout.write("\n\n\t> Converting plink data to binary format for " + sample_id)
                print("plink --file " + sample_id + " --out " + sample_id + " --make-bed --allow-no-sex")
                os.system("plink --file " + sample_id + " --out " + sample_id + " --make-bed --allow-no-sex")
                fout.write(time.strftime(
                    "\n\t%H:%M:%S %Z ") + "\tplink --file " + sample_id + " --out " + sample_id + " --make-bed --allow-no-sex")
                ## Move new binary files to "bin4pca" folder
                for item in os.listdir(cwd):
                    if sample_id + ".bed" in item:
                        os.system('mv ' + item + ' ' + cwd + '/' + "bin4pca")
                for item in os.listdir(cwd):
                    if sample_id + ".bim" in item:
                        os.system('mv ' + item + ' ' + cwd + '/' + "bin4pca")
                for item in os.listdir(cwd):
                    if sample_id + ".fam" in item:
                        os.system('mv ' + item + ' ' + cwd + '/' + "bin4pca")
                        # for item in os.listdir(cwd):
                        #     if sample_id+".log" in item:
                        #         os.system('mv '+item+' '+cwd+'/'+"bin4pca")
    elif conv2bin == False:
        pass

    ## Merge dataset and our sample
    if mergeIt == True:
        ## Do reduce dataset's SNP list to match the sample's
        if reduceDSsnp == True:
            cwdBIN = cwd + "/bin4pca"
            os.chdir(cwdBIN)
            datasetName = mapFile.strip(".bim")
            datasetPath = mapFileDatasetPath.strip(".bim")
            for i in os.listdir(cwdBIN):
                if ".bed" == i[-4:] and datasetName not in i and "excluded" not in i and "mergedex" not in i:  ## SAMPLE NAME CANNOT INCLUDE "EXCLUDED"
                    sample_id, sep, rest = i.partition('.bed')
                    print("\t> Creating SNP list from sample " + sample_id)
                    fout.write("\n\n\t> Creating SNP list from sample " + sample_id)

                    script = open("reduceSNPlist.R", 'w')
                    script.write("## Create list of SNPs from BIM file\n")
                    script.write('setwd("')
                    script.write(cwdBIN)
                    script.write('")\nwdfiles=list.files(getwd())')
                    script.write('\n\nfor(file in wdfiles) {\n\tif(grepl(".bim$",file)==TRUE & grepl("excluded",file)==FALSE & grepl("reduced",file)==FALSE & grepl("mergedex",file)==FALSE) {')
                    script.write('\n\t\t#print(file)\n\t\tsample_id=strsplit(file,".bim")\n\t\t#print(sample_id)\n\t\tfile_read=read.table(as.character(file), as.is=1)')
                    script.write('\n\t\tsnpList=file_read$V2\n\t\twrite.table(snpList, paste(sample_id,"_snpList.txt",sep=""), quote = FALSE,col.names = FALSE, row.names = FALSE)\n\t}\n}')
                    script.close()
                    os.system("R CMD BATCH reduceSNPlist.R")
                    fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tR CMD BATCH reduceSNPlist.R")

                    print("\t> Reducing dataset's SNPs to match the sample's")
                    print("plink --bfile " + datasetPath + " --extract " + sample_id + "_snpList.txt --make-bed --out " + datasetName + "_" + sample_id + "_reduced --allow-no-sex")  ## dataPP_D1_reduced.bed/bim/fam
                    os.system("plink --bfile " + datasetPath + " --extract " + sample_id + "_snpList.txt --make-bed --out " + datasetName + "_" + sample_id + "_reduced --allow-no-sex")
                    fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tplink --bfile " + datasetPath + " --extract " + sample_id + "_snpList.txt --make-bed --out " + datasetName + "_" + sample_id + "_reduced --allow-no-sex")

                    print("\t\t\tplink --bfile " + datasetName + "_" + sample_id + "_reduced --bmerge " + sample_id + " --make-bed --out " + sample_id + " --allow-no-sex")  ## D1-merge.fam + D1-merge.missnp
                    os.system("plink --bfile " + datasetName + "_" + sample_id + "_reduced --bmerge " + sample_id + " --make-bed --out " + sample_id + " --allow-no-sex")
                    fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tplink --bfile " + datasetName + "_" + sample_id + "_reduced --bmerge " + sample_id + " --make-bed --out " + sample_id + " --allow-no-sex")

                    if sample_id + "-merge.missnp" in os.listdir(cwdBIN):
                        print("\t > Found missing SNPsfile, excluding...")
                        print("\t\t\t\tplink --bfile " + sample_id + " --exclude " + sample_id + "-merge.missnp --make-bed --out " + sample_id + "excluded --allow-no-sex")  ## D1excluded.bed/bim/fam
                        os.system("plink --bfile " + sample_id + " --exclude " + sample_id + "-merge.missnp --make-bed --out " + sample_id + "excluded --allow-no-sex")
                        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tplink --bfile " + sample_id + " --exclude " + sample_id + "-merge.missnp --make-bed --out " + sample_id + "excluded --allow-no-sex")

                        print("\t\t\t\tplink --bfile " + datasetName + "_" + sample_id + "_reduced --exclude " + sample_id + "-merge.missnp --make-bed --out " + datasetName + "red_" + sample_id + "excluded --allow-no-sex")  ## dataPPred_D1excluded.bed/bim/fam
                        os.system("plink --bfile " + datasetName + "_" + sample_id + "_reduced --exclude " + sample_id + "-merge.missnp --make-bed --out " + datasetName + "red_" + sample_id + "excluded --allow-no-sex")
                        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tplink --bfile " + datasetName + "_" + sample_id + "_reduced --exclude " + sample_id + "-merge.missnp --make-bed --out " + datasetName + "red_" + sample_id + "excluded --allow-no-sex")

                        print("\t\t\t\tplink --bfile " + datasetName + "red_" + sample_id + "excluded --bmerge " + sample_id + "excluded" + " --make-bed --out " + sample_id + "mergedex --allow-no-sex")  ## D1mergedex.bed/bim/fam
                        os.system("plink --bfile " + datasetName + "red_" + sample_id + "excluded --bmerge " + sample_id + "excluded" + " --make-bed --out " + sample_id + "mergedex --allow-no-sex")
                        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tplink --bfile " + datasetName + "red_" + sample_id + "excluded --bmerge " + sample_id + "excluded" + " --make-bed --out " + sample_id + "mergedex --allow-no-sex")

                        ## Apply missing genotype rate SNP exclusion step
                        if genoPlink != 1:
                            print("\t\t\t\tplink --bfile " + sample_id + "mergedex --geno " + str(genoPlink) + " --make-bed --out " + sample_id + "mergedex --allow-no-sex")  ## D1mergedex.bed/bim/fam
                            os.system("plink --bfile " + sample_id + "mergedex --geno " + str(genoPlink) + " --make-bed --out " + sample_id + "mergedex --allow-no-sex")
                            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tplink --bfile " + sample_id + "mergedex --geno " + str(genoPlink) + " --make-bed --out " + sample_id + "mergedex --allow-no-sex")

                        for i in os.listdir(cwdBIN):
                            if sample_id + "mergedex.bed" == i:
                                sample_id, sep, rest = i.partition('.bed')
                                print(sample_id)
                                convPed2Eig = open(sample_id + "_CONV_PED_2_EIG", 'w')
                                convPed2Eig.write("genotypename:\t" + sample_id + ".bed\nsnpname:\t" + sample_id + ".bim\nindivname:\t" + sample_id + ".fam\noutputformat:\tEIGENSTRAT\ngenotypeoutname:\t" + sample_id + ".eigenstratgeno\nsnpoutname:\t" + sample_id + ".snp\nindivoutname:\t" + sample_id + ".ind\nfamilyname:\tYES")
                                convPed2Eig.close()
                                print("convertf -p " + sample_id + "_CONV_PED_2_EIG")
                                os.system("convertf -p " + sample_id + "_CONV_PED_2_EIG")
                                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tconvertf -p " + sample_id + "_CONV_PED_2_EIG")

                    else:  ## Make sure files in this case are outputed as XXX-merge.evec
                        print(sample_id)
                        sample_id_complete = sample_id + "-merge"
                        print(sample_id_complete)
                        os.rename(str(sample_id + ".fam"), str(sample_id + "-merge.fam"))
                        os.rename(str(sample_id + ".bim"), str(sample_id + "-merge.bim"))
                        os.rename(str(sample_id + ".bed"), str(sample_id + "-merge.bed"))
                        convPed2Eig = open(sample_id_complete + "_CONV_PED_2_EIG", 'w')
                        convPed2Eig.write("genotypename:\t" + sample_id_complete + ".bed\nsnpname:\t" + sample_id_complete + ".bim\nindivname:\t" + sample_id_complete + ".fam\noutputformat:\tEIGENSTRAT\ngenotypeoutname:\t" + sample_id_complete + ".eigenstratgeno\nsnpoutname:\t" + sample_id_complete + ".snp\nindivoutname:\t" + sample_id_complete + ".ind\nfamilyname:\tYES")
                        convPed2Eig.close()
                        print("convertf -p " + sample_id_complete + "-merge_CONV_PED_2_EIG")
                        os.system("convertf -p " + sample_id + "-merge_CONV_PED_2_EIG")
                        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tconvertf -p " + sample_id + "-merge_CONV_PED_2_EIG")

                        # Now go to R and convert the mergedex.IND to be the same as the original IND file in terms of columns and populations

            print("\t> Correcting format of *.ind files")
            fout.write("\n\n\t> Correcting format of *.ind files")
            script = open("pca_indfile_correction.R", 'w')
            script.write('setwd("')
            script.write(cwdBIN)
            script.write('")\nwdfiles=list.files(getwd())\n\nfor(i in wdfiles) {\n\tif(grepl(".ind$",i)==TRUE) {')
            script.write('\n\t\tname=strsplit(i, ".ind")\n\t\tmerged_ind=read.table(i, as.is=1, stringsAsFactors = FALSE, sep=":")\n\t\tmerged_ind_temp1=data.frame(as.character(merged_ind[,2]))')
            script.write('\n\t\twrite.table(merged_ind_temp1, "merged_ind_temp1.txt", sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE)\n\t\tmerged_ind_temp2=read.table("merged_ind_temp1.txt", sep=" ")')
            script.write('\n\t\tmerged_ind_temp2$V3={}\n\t\tmerged_ind_temp2$V4={}\n\t\tmerged_ind_temp2$V5={}\n\t\tmerged_ind_temp2$V6=merged_ind$V1\n\t\tfinal_name=as.character(i)\n\t\twrite.table(merged_ind_temp2, final_name, sep=" ",quote=FALSE,col.names=FALSE,row.names=FALSE)')
            script.write('\n\t}\n}')
            script.close()
            os.system("R CMD BATCH pca_indfile_correction.R")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tR CMD BATCH pca_indfile_correction.R")

    elif mergeIt == False:
        pass

    ## Create PCA parameter files to be used by smartpca
    if PCAparamFile == True:
        cwdBIN = cwd + "/bin4pca"
        os.chdir(cwdBIN)
        listForPCA = []
        print("\t> Creating PCA parameter files to be used by smartpca")
        fout.write("\n\n\t> Creating PCA parameter files to be used by smartpca")
        for i in os.listdir(cwdBIN):
            if ".eigenstratgeno" in i:
                listForPCA.append(i)
                sample_id, sep, rest = i.partition('.eigenstratgeno')
                createPCA = open(sample_id + "_CREATE_PCA", 'w')
                if seqPlatform == "MiSeq":
                    createPCA.write(
                        "genotypename:\t" + i + "\nsnpname:\t" + sample_id + ".snp\nindivname:\t" + sample_id + ".ind\nevecoutname:\t" + sample_id + ".evec\nevaloutname:\t" + sample_id + ".eval")
                    if popList != None:
                        createPCA.write("\npoplistname:\t" + popList)
                    createPCA.write("\nnumoutevec:\t4\nnumoutlieriter:\t0\nkillr2:\tYES\nr2thresh:\t0.2")

                    ## Not sure about this step here. It might not be needed to project when reducing all SNPs
                    ## Latest info says it should always be projected, even for low coverage 1-2X or 340K SNP capture
                    if reduceDSsnp == True:
                        createPCA.write("\nlsqproject:\tYES")
                    if reduceDSsnp == False:
                        createPCA.write("\nlsqproject:\tYES")
                    createPCA.close()
                    print("smartpca -p " + sample_id + "_CREATE_PCA")
                    os.system("smartpca -p " + sample_id + "_CREATE_PCA")
                    fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tsmartpca -p " + sample_id + "_CREATE_PCA")
                elif seqPlatform == "NextSeq":
                    createPCA.write(
                        "genotypename:\t" + i + "\nsnpname:\t" + sample_id + ".snp\nindivname:\t" + sample_id + ".ind\nevecoutname:\t" + sample_id + ".evec\nevaloutname:\t" + sample_id + ".eval")
                    if popList != None:
                        createPCA.write("\npoplistname:\t" + popList)
                    createPCA.write("\nnumoutevec:\t4\nnumoutlieriter:\t0\nr2thresh:\t0.2")

                    ## Not sure about this step here. It might not be needed to project when reducing all SNPs
                    ## Latest info says it should always be projected, even for low coverage 1-2X or 340K SNP capture
                    if reduceDSsnp == True:
                        createPCA.write("\nlsqproject:\tYES")
                    if reduceDSsnp == False:
                        createPCA.write("\nlsqproject:\tYES")
                    createPCA.close()
                    print("smartpca -p " + sample_id + "_CREATE_PCA")
                    os.system("smartpca -p " + sample_id + "_CREATE_PCA")
                    fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tsmartpca -p " + sample_id + "_CREATE_PCA")

    if plot != False:  ## Might give UnboundLocalError if PCAparamfile=False
        print("\t> Plotting low-res PCA")
        fout.write("\n\n\t> Plotting low-res PCA")
        script = open("PCA_plot.R", 'w')
        script.write('setwd("')
        script.write(cwdBIN)
        script.write('")\nwdfiles=list.files(getwd()')
        script.write(')\n\npop_db=read.csv("')
        script.write(plot)
        script.write('", sep=",", header=FALSE, stringsAsFactors=FALSE)\n')
        script.write('for(file in wdfiles) {\n\tif(grepl(".evec$",file)==TRUE) {\n\t\tsample_id=strsplit(file,".evec")\n\t\tfile_read=read.table(as.character(file), as.is=1,row.names = NULL)\n\t\tif(grepl("mergedex",file)==TRUE) {')
        script.write('\n\t\t\ttrue_id=strsplit(file,"mergedex.evec")\n\t\t\tsnp_count=read.table(paste(as.character(true_id),"mergedex.snp", sep=""))\n\t\t\tsnp_count=length(snp_count$V1)\n\t\t}')
        script.write('\n\t\tif(grepl("-merge",file)==TRUE) {\n\t\t\ttrue_id=strsplit(file,"-merge.evec")\n\t\t\tsnp_count=read.table(paste(as.character(true_id),"-merge.snp", sep=""))\n\t\t\tsnp_count=length(snp_count$V1)\n\t\t}')
        script.write('\n\t\tpdf(file=paste(getwd(),"/",as.character(true_id),".pdf",sep=""),pointsize=30,width=40,height=30)')
        script.write('\n\t\tplot(1,type="n",ylab="PC2",xlab="PC1",ylim=c(min(file_read[,3]),max(file_read[,3])),xlim=c(min(file_read[,2]),max(file_read[,2])))')
        script.write('\n\t\trows=as.numeric(length(file_read$V1))')
        script.write('\n\t\tit=0')
        script.write('\n\t\tfor(row in rows) {')
        script.write('\n\t\t\twhile(it<length(file_read$V1)) {\n\t\t\t\tit=it+1\n\t\t\t\tif(file_read[it,1] %in% pop_db$V2==TRUE) {')
        script.write('\n\t\t\t\t\tpos_in_pop_db=which(pop_db$V2==file_read[it,1])\n\t\t\t\t\ttext(file_read[it,2], file_read[it,3], pop_db$V10[pos_in_pop_db], col=rgb(120,120,120,100,maxColorValue=255),cex=0.7,font=2)')
        script.write('\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\n\t\tpoints(file_read[which(file_read$V1==true_id),2],file_read[which(file_read$V1==true_id),3], bg="firebrick2", cex=1.5,pch = 21)')
        script.write('\n\t\tlegend("bottomleft", legend=c(paste0(true_id," : ",snp_count," SNPs")), text.col="black", pch=21,pt.bg="firebrick2", y.intersp = 0.9,pt.cex = 2,bty="n")\n\t\tdev.off()\n\t}\n}')
        script.close()
        os.system("R CMD BATCH PCA_plot.R")
        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tR CMD BATCH PCA_plot.R")
        os.remove(cwdBIN+"/PCA_plot.R")
        os.remove(cwdBIN + "/PCA_plot.Rout")


        # TODO
        ## In final PCA report explain what kind of PCA was plotted according to the cmd line arguments of the function

        ## PED
        # Family ID
        # Sample ID
        # Paternal ID
        # Maternal ID
        # Sex (1=male; 2=female; other=unknown)
        # Affection (0=unknown; 1=unaffected; 2=affected)
        # Genotypes (space or tab separated, 2 for each marker. 0=missing)

        # FamilyID  SampleID   PaternalID   MaternalID  Sex    Affection    Genotypes
        # FAM1      NA06985         0           0        1          1           A     T     T     T     G     G     C     C     A     T     T     T     G     G     C     C
        # FAM1      NA06991         0           0        1          1           C     T     T     T     G     G     C     C     C     T     T     T     G     G     C     C

        ## MAP
        # 21     rs11511647     0     26765
        # X     rs3883674     0     32380



        ## SAMTOOLS PILEUP FORMAT DESCRIPTION ##
        #    chr2           41784              A            1               ,$                   A
        # chromosome name, coordinate, reference base, read bases, read qualities and alignment mapping qualities
        #   col1            col2             col3       col4            col5                col6

        ## col5 Read qualities INFO ##
        #  .  - match on reference base on + strand
        #  ,  - match on reference base on - strand
        #  >/<  - reference skip
        #  ACTGN  - mismatch on + strand
        #  actgn  - mismatch on - strand
        #  \+[0-9]+[ACGTNacgtn]+  - insertion of length [0-9] between this ref position and next
        #  -[0-9]+[ACGTNacgtn]+  - deletion from ref. Bases presented as * in following lines
        #  ^  - start of read. ASCII of the character following ?^? minus 33 gives the mapping quality.
        #  $  - end of read

        ## col6 Alihnment mapping qualities INFO ##
        #  A  - character in ASCII format with corresponding phred quality (ASCII-33)


        ## GATK PILEUP FORMAT DESCRIPTION ##
        #    chr2           41784              A                AAAAATA                                    F
        # chromosome name, coordinate, reference base, read bases (7 reads in example above), and base quality
        #   col1            col2             col3           col4                                   col5

        ## col5 Alihnment mapping qualities INFO ##


#  A  - character in ASCII format with corresponding phred quality (ASCII-33)


def moveFilesBatch(byExtension, bySample=False, copy=False, singleFolder=False, asRegex=False):
    """
    Moves files into new folder(s). Works by extension of the files, and/or by sample identifier at the start of the files. Please note that input is as regular expressions and include special characters.
    NO REGEX AFTER ALL - TRUE! DO NOT USE REGEX FOR INPUTTING FILE EXTENSIONS OR SAMPLE NAMES

    ############### WRONG SECTION
    For information on how to use special characters, please check https://docs.python.org/3.5/library/re.html

    Note on Special Characters:
        All input must be in the form of regular expressions, accepting special characters. This allows for more accurate tasks, as in the two examples below.

        # Using the special character '.':
        The dot in a regular expression matches any character, except a newline '\n'.
        >>> moveFilesBatch(byExtension = False, bySample = "foo_L00.", copy = True, singleFolder=False)
        This will retrieve all files named "foo_L00.", with any character matching the dot, as in "foo_L001" and "foo_L002".

        # Using the special character '\':
        The backslash escapes special characters (allowing to match characters like '*', '.', and so forth).
        >>> moveFilesBatch(byExtension = "^foo\.bar\.tar", bySample = False, copy = True, singleFolder=False)
        This will retrieve all files starting with "foo.bar.tar". If backslash is not used, all files starting with "fooXbarXtar" would be retrieved, as '.' would be considered the special character explained above.

    Note on Multiple Inputs:
        To input more than one extension/sample please use separator '/', as below:
        >>> moveFilesBatch(byExtension = "\.foo", bySample = False, copy = True, singleFolder=False)
        >>> moveFilesBatch(byExtension = "\.foo/\.bar", bySample = False, copy = True, singleFolder=False)
    ################

    :param byExtension:
    (boolean-False)
    (string)    Boolean or string with desired extension to match at the end of the files. True boolean is not accepted

    :param bySample:
    (boolean-False)
    (string)    Boolean or string with desired text to match at the start of files. True boolean is not accepted - a sample name (or list of names) must be specified
    This parameter is exclusive to use with NextSeq data and files should be in the form "foo_Sxx_Lxxx*" in order to identify each sample correctly.

    :param copy:
    (boolean)   Whether to copy (True) or move files (False)

    :param singleFolder:
    (boolean-False)
    (string)    Boolean or string with desired text for folder name. True bollean is not accepted - a folder name must be specified

    :param asRegex:
    (boolean)   Whether the input strings are in regular expression format or not. False by default, matching string literals instead

    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< BATCHMOVING - Moving specific files to new folder(s) >")
    fout.write("\n< BATCHMOVING - Moving specific files to new folder(s) >")

    ## Argument checking
    if byExtension == False and bySample == False and singleFolder == False:
        print("\t> ERROR: Please check your arguments. No input was provided, all arguments set as False.")
        exit()
    if byExtension != False or bySample != False or singleFolder != False:
        if byExtension == True or bySample == True or singleFolder == True:
            print(
                "\t> ERROR: Please check your arguments. 'byExtension', 'bySample', and 'singleFolder' must not be given the boolean True.")
            exit()
        paddy1 = 0
        if type(byExtension) != str and byExtension != False:
            paddy1 += 1
        if type(bySample) != str and bySample != False:
            paddy1 += 1
        if type(singleFolder) != str and singleFolder != False:
            paddy1 += 1
        if paddy1 != 0:
            print(
                "\t> ERROR: Please check your arguments. 'byExtension', 'bySample', and 'singleFolder' must either be False or a string.")
            exit()

    ## Checking if multiple input exists for byExtension and/or bySample, recognized by separator '/'
    if type(byExtension) == str and "/" in byExtension:
        # byExtension = byExtension.split("/")
        print("A =", byExtension)
        byExtension = re.split("/", byExtension)
        print("B =", byExtension)
        # byExtension[0] = byExtension[0].replace(r"\\\\","")
        print("C =", byExtension)
        print("ByExt0:", repr(byExtension))

    elif byExtension != False:
        # print(byExtension)
        byExtension = [byExtension]
        print("ByExt1:", byExtension)

    if type(bySample) == str and "/" in bySample:
        bySample = bySample.split("/")
        print("BySample1:", bySample)
    elif bySample != False:
        bySample = [bySample]
        print(bySample)

        # print(re.compile(byExtension[0]))

        # a = "abc.abc.abc"
        # b = r"abc.abc.abc"
        # c = re.escape("abc.abc.abc")
        # print(a)
        # print(b)
        # print(c)
        # print(re.escape("abc.abc/abc"))

    ## For bySample = False and byExtension = string (single or multi)
    ## Search for files with desired extensions and add to dictionary with key as the file and value as the extension
    if bySample == False and byExtension != False:
        fileExtDic = {}
        for i in os.listdir(cwd):
            if len(byExtension) == 1:
                if re.search(byExtension[0] + "$", i) != None and os.path.isdir(i) != True:
                    fileExtDic[i] = byExtension[0]
            elif len(byExtension) > 1:
                for q in byExtension:
                    if re.search(q + "$", i) != None and os.path.isdir(i) != True:
                        print(q)
                        q = q.replace("\.", ".")
                        print(q)
                        fileExtDic[
                            i] = q  ###### ???? Maybe not anymore #### PROBLEM HERE, NEED TO FIND A WAY TO TAKE THE '\' AND ANY CONSEadaIVE CHARACTER FROM STRING. SHOULD BE OK USING THE INDEX FOR THE POSITION OF THE '\' AND THEN +1

        ## Creating new folder(s)
        if type(singleFolder) == str:  ## Creating single folder
            newdir = cwd + '/' + singleFolder
            if singleFolder in os.listdir(cwd):
                pass
            else:
                os.mkdir(newdir)
        elif singleFolder == False:  ## Creating one folder per extension
            if len(byExtension) == 1:
                newdir = cwd + '/' + byExtension[0]
                if byExtension[0] in os.listdir(cwd):
                    pass
                else:
                    os.mkdir(newdir)
            elif len(byExtension) > 1:
                for q in byExtension:
                    newdir = cwd + '/' + q
                    if q in os.listdir(cwd):
                        pass
                    else:
                        os.mkdir(newdir)

        print(fileExtDic)
        ## Moving or copying files
        ## Move all to one single folder specified in singleFolder argument
        if singleFolder != False and type(singleFolder) is str:
            newdir = cwd + '/' + singleFolder
            for i in fileExtDic:
                if copy == False:  ## Move files
                    print('mv -b -f ' + cwd + '/' + i + ' ' + newdir + '/')
                    os.system('mv -b -f ' + cwd + '/' + i + ' ' + newdir + '/')
                else:  ## Copy files
                    print('cp -b -f ' + cwd + '/' + i + ' ' + newdir + '/')
                    os.system('cp -b -f ' + cwd + '/' + i + ' ' + newdir + '/')
        else:  ## Move to folder specified in byExtension
            for i in fileExtDic:
                if copy == False:  ## Move files
                    print('mv -b -f ' + cwd + '/' + i + ' ' + cwd + '/' + fileExtDic[i] + '/')
                    os.system('mv -b -f ' + cwd + '/' + i + ' ' + cwd + '/' + fileExtDic[i] + '/')
                else:  ## Copy files
                    print('cp -b -f ' + cwd + '/' + i + ' ' + cwd + '/' + fileExtDic[i] + '/')
                    os.system('cp -b -f ' + cwd + '/' + i + ' ' + cwd + '/' + fileExtDic[i] + '/')



                    ## For bySample = string (single or multi) and byExtension = False
                    ## For automatically moving all files from each sample to individual folders
    ## Searches for files from the desired sample(s) and add to dictionary with key as the file and value as the sample name
    if bySample != False and byExtension == False:
        sampleDic = {}
        for i in os.listdir(cwd):
            if len(bySample) == 1:
                if re.search("^" + bySample[0], i) != None and os.path.isdir(i) != True:
                    sampleDic[i] = bySample[0]
            elif len(bySample) > 1:
                for q in bySample:
                    stri = "^" + q + "(?=_S\d\d_L\d\d)"
                    if re.search(stri, i) != None and os.path.isdir(i) != True:
                        sampleDic[i] = q

        ## Creating new folder(s)
        if type(singleFolder) == str:  ## Creating single folder
            newdir = cwd + '/' + singleFolder
            if singleFolder in os.listdir(cwd):
                pass
            else:
                os.mkdir(newdir)
        elif singleFolder == False:  ## Creating one folder per sample
            if len(bySample) == 1:
                newdir = cwd + '/' + bySample[0]
                if bySample[0] in os.listdir(cwd):
                    pass
                else:
                    os.mkdir(newdir)
            elif len(bySample) > 1:
                for q in bySample:
                    newdir = cwd + '/' + q
                    if q in os.listdir(cwd):
                        pass
                    else:
                        os.mkdir(newdir)

                        ## Moving or copying files
                        ## Move all to one single folder specified in singleFolder argument
        if singleFolder != False and type(singleFolder) is str:
            newdir = cwd + '/' + singleFolder
            for i in sampleDic:
                if copy == False:  ## Move files
                    print('mv -b -f ' + cwd + '/' + i + ' ' + newdir + '/')
                    # os.system('mv -b -f '+cwd+'/'+i+' '+newdir+'/')
                else:  ## Copy files
                    print('cp -b -f ' + cwd + '/' + i + ' ' + newdir + '/')
                    # os.system('cp -b -f '+cwd+'/'+i+' '+newdir+'/')
        else:  ## Move to individual folder specified in byExtension
            for i in sampleDic:
                if copy == False:  ## Move files
                    print('mv -b -f ' + cwd + '/' + i + ' ' + cwd + '/' + sampleDic[i] + '/')
                    # os.system('mv -b -f '+cwd+'/'+i+' '+cwd+'/'+sampleDic[i]+'/')
                else:  ## Copy files
                    print('cp -b -f ' + cwd + '/' + i + ' ' + cwd + '/' + sampleDic[i] + '/')
                    # os.system('cp -b -f '+cwd+'/'+i+' '+cwd+'/'+sampleDic[i]+'/')

        print(sampleDic)

    ## For bySample and byExtension
    if bySample != False and byExtension != False:
        sampleDic = {}
        for i in os.listdir(cwd):
            if len(bySample) == 1:
                if re.search("^" + bySample[0], i) != None and os.path.isdir(i) != True:
                    sampleDic[i] = bySample[0]
            elif len(bySample) > 1:
                for q in bySample:
                    stri = "^" + q + "(?=_S\d\d_L\d\d)"
                    if re.search(stri, i) != None and os.path.isdir(i) != True:
                        print(i)
                        for o in byExtension:
                            strin = o + "$"
                            print(strin)
                            print(i)
                            if re.search(strin, i) != None:
                                sampleDic[i] = q

        ## Creating new folder(s)
        if type(singleFolder) == str:  ## Creating single folder
            newdir = cwd + '/' + singleFolder
            if singleFolder in os.listdir(cwd):
                pass
            else:
                os.mkdir(newdir)
        elif singleFolder == False:  ## Creating one folder per sample
            if len(bySample) == 1:
                newdir = cwd + '/' + bySample[0]
                if bySample[0] in os.listdir(cwd):
                    pass
                else:
                    os.mkdir(newdir)
            elif len(bySample) > 1:
                for q in bySample:
                    newdir = cwd + '/' + q
                    if q in os.listdir(cwd):
                        pass
                    else:
                        os.mkdir(newdir)

                        ## Moving or copying files
                        ## Move all to one single folder specified in singleFolder argument
        if singleFolder != False and type(singleFolder) is str:
            newdir = cwd + '/' + singleFolder
            for i in sampleDic:
                if copy == False:  ## Move files
                    print('mv -b -f ' + cwd + '/' + i + ' ' + newdir + '/')
                    os.system('mv -b -f ' + cwd + '/' + i + ' ' + newdir + '/')
                else:  ## Copy files
                    print('cp -b -f ' + cwd + '/' + i + ' ' + newdir + '/')
                    os.system('cp -b -f ' + cwd + '/' + i + ' ' + newdir + '/')
        else:  ## Move to individual folder specified in byExtension
            for i in sampleDic:
                if copy == False:  ## Move files
                    print('mv -b -f ' + cwd + '/' + i + ' ' + cwd + '/' + sampleDic[i] + '/')
                    os.system('mv -b -f ' + cwd + '/' + i + ' ' + cwd + '/' + sampleDic[i] + '/')
                else:  ## Copy files
                    print('cp -b -f ' + cwd + '/' + i + ' ' + cwd + '/' + sampleDic[i] + '/')
                    os.system('cp -b -f ' + cwd + '/' + i + ' ' + cwd + '/' + sampleDic[i] + '/')

        print(sampleDic)


def mergeBams(sample_id="Output", qualThreshold = 30, version = "1.3.1", samtoolsOrGATK="samtools"):
    """
    Merges **ALL** rmdup bam files in folder.

    Input:
    Foo1.q30.rmdup.bam
    Foo2.q30.rmdup.bam
    Foo3.q30.rmdup.bam

    Output:
    Foo_merged.q30.rmdup.bam
    Foo_merged.q30.rmdup.bai

    :param sample_id:
    (string)    Identifier for output file. "Output" if none provided.

    :param qualThreshold:
    (int)   Quality score (QC) filter for the newly merged files. Default = 30.

    :param version:
    (str)   Samtools' version, either 0.1.x or 1.x.x.

    :param samtoolsOrGATK:
    (string)    Use SAMTOOLS or GATK to merge the files. Accepts "samtools" and "gatk"/"GATK".
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< MERGEBAMS - Merging all bam files in current folder >")
    fout.write("\n< MERGEBAMS - Merging all bam files in current folder >")

    bamList = []
    for i in os.listdir(cwd):
        if '.rmdup.bam' == i[-10:]:
            bamList.append(i)

    qualThreshold = str(qualThreshold)

    command = "picard.jar MergeSamFiles "
    if len(bamList) != 0:
        print("\t> Merging all bam files for sample with id " + sample_id + " and sorting")
        fout.write("\n\n\t> Merging all bam files for sample with id " + sample_id + " and sorting")
        for i in bamList:
            command = command + "I=" + i + " "
    command = command + "O=" + sample_id + "_merged.q00.bam USE_THREADING=True"
    print(command)
    os.system(command)
    fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + command)

    for i in os.listdir(cwd):
        if "_merged.q00.bam" == i[-15:]:
            print("\tsamtools index " + i)
            os.system("samtools index " + i)

            if version[0:3] == "0.1":
                sample_id_loc = i[:-4]
                print("\t" + 'samtools sort ' + i + ' ' + sample_id_loc + '.sort' + "\n")
                os.system('samtools sort ' + i + ' ' + sample_id_loc + '.sort')

            if version[0:2] == "1.":
                sample_id_loc = i[:-4]
                print("\t" + 'samtools sort ' + i + ' -o ' + sample_id_loc + '.sort.bam' + "\n")
                os.system('samtools sort ' + i + ' -o ' + sample_id_loc + '.sort.bam')

            newI = sample_id_loc + '.sort.bam'

            print("\tsamtools index " + newI)
            os.system("samtools index " + newI)

            print("\tsamtools view -q " + qualThreshold + " -b " + newI + " | samtools rmdup -s " + newI + " " + sample_id + "_merged.q" + qualThreshold + ".rmdup.bam")
            os.system("samtools view -q " + qualThreshold + " -b " + newI + " | samtools rmdup -s " + newI + " " + sample_id + "_merged.q" + qualThreshold + ".rmdup.bam")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools view -q " + qualThreshold + " -b " + newI + " | samtools rmdup -s " + newI + " " + sample_id + "_merged.q" + qualThreshold + ".rmdup.bam")

            newI2_full = sample_id + "_merged.q" + qualThreshold + ".rmdup.bam"
            newI2 = sample_id + "_merged.q" + qualThreshold
            newI3 = sample_id + "_merged"

            print("\tsamtools flagstat " + newI2_full + " > " + newI2 + ".summary.txt")
            os.system("samtools flagstat " + newI2_full + " > " + newI2 + ".summary.txt")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools flagstat " + newI2_full + " > " + newI2 + ".summary.txt")
            print("\tsamtools index " + newI2_full + " " + newI2_full + ".bai")
            os.system("samtools index " + newI2_full + " > " + newI2_full + ".bai")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools index " + newI2_full + " > " + newI2_full + ".bai")
            print("\tsamtools idxstats " + newI2_full + " > " + newI3 + ".chr.txt")
            os.system("samtools idxstats " + newI2_full + " > " + newI3 + ".chr.txt")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools idxstats " + newI2_full + " > " + newI3 + ".chr.txt")


def coverage(contigFileLocation):
    """
    Calculates coverage for all *.bam files with removed duplicates in folder, returning the values in the terminal and in "Coverages.txt".
    Requires contig size file with the information in the following format:
    chr1    249250621
    chr2    243199373
    (...)

    Input:
    Foo.q30.rmdup.bam
    or
    Bar_merged.q30.rmdup.bam

    Output:
    Foo_coverage.txt

    :param contigFileLocation:
    (string)   Path to text file containing contig size information
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< COVERAGE - Assess genome coverage >")
    fout.write("\n< COVERAGE - Assess genome coverage >")

    for i in os.listdir(cwd):
        if "_merged" in i:
            if re.search("_merged\.q..\.rmdup\.bam$", i) != None:
            #if "_merged.q30.rmdup.bam" == i[-21:]:
                sample_id, sep, rest = i.partition("_merged.")
                print("\t> Creating file with coverage depths for " + sample_id)
                fout.write("\n\n\t> Creating file with coverage depths for " + sample_id)
                print(
                    "\tgenomeCoverageBed -ibam " + i + " -g " + contigFileLocation + " > " + sample_id + "_coverage.txt")
                os.system(
                    "genomeCoverageBed -ibam " + i + " -g " + contigFileLocation + " > " + sample_id + "_coverage.txt")
                fout.write(time.strftime(
                    "\n\t%H:%M:%S %Z ") + "\t" + "genomeCoverageBed -ibam " + i + " -g " + contigFileLocation + " > " + sample_id + "_coverage.txt")
            else:
                pass
        elif re.search("\.q..\.rmdup\.bam$", i) != None and "merge" not in i:
        #elif ".q30.rmdup.bam" == i[-14:] and "merge" not in i:
            #sample_id, sep, rest = i.partition(".q30")
            sample_id, rest = re.split('\.q..', i)
            sep = re.findall('\.q..', i)[0]

            print("\t> Creating file with coverage depths for " + sample_id)
            fout.write("\n\n\t> Creating file with coverage depths for " + sample_id)
            print("\tgenomeCoverageBed -ibam " + i + " -g " + contigFileLocation + " > " + sample_id + "_coverage.txt")
            os.system(
                "genomeCoverageBed -ibam " + i + " -g " + contigFileLocation + " > " + sample_id + "_coverage.txt")
            fout.write(time.strftime(
                "\n\t%H:%M:%S %Z ") + "\t" + "genomeCoverageBed -ibam " + i + " -g " + contigFileLocation + " > " + sample_id + "_coverage.txt")

    for i in os.listdir(cwd):
        if "_coverage.txt" in i:
            coverages = open("Coverages.txt", "w")
            coverages.write("Sample\tGenome Covered\tGenome Depth\tMT covered\tMT Depth\n")

    for i in os.listdir(cwd):
        ## Check coverage for whole genome
        if "_coverage.txt" in i:
            sample_id, sep, rest = i.partition("_coverage")
            print("\n\t> Calculating coverage depths for " + sample_id + " from " + i)
            fout.write("\n\n\t> Calculating coverage depths for " + sample_id + " from " + i)
            coverFile = open(i, "r")
            prod1 = 0
            sum2 = 0
            for line in coverFile:
                if "genome" in line:
                    chromo, depth, numb, length, perc = line.split("\t")
                    if depth == "0":
                        genomeCovered = "%.2f" % float((1 - float(perc)) * 100)
                    prod1 = prod1 + (float(depth) * float(numb))
                    sum2 = sum2 + float(numb)
            genCovDepthFinal = float(prod1) / float(sum2)
            genCovDepthFinal = "{0:.4f}".format(float(genCovDepthFinal))
            if float(genCovDepthFinal) < 0.0001:
                genCovDepthFinal = "Below 0.0001"
            print("\tGenome covered: " + genomeCovered + "%")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "Genome covered: " + genomeCovered + "%")
            print("\tGenome coverage depth: " + genCovDepthFinal + "X")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "Genome coverage depth: " + genCovDepthFinal + "X")

            ## Check coverage for mitochondrial genome
            coverFile.seek(0)
            prodMT = 0
            sumMT = 0
            mtCovered = 0
            for line in coverFile:
                if "chrM" in line:
                    chromo, depth, numb, length, perc = line.split("\t")
                    if depth == "0":
                        mtCovered = "%.2f" % float((1 - float(perc)) * 100)
                    prodMT = prodMT + (float(depth) * float(numb))
                    sumMT = sumMT + float(numb)

            if sumMT >0:
                mtCovDepthFinal = float(prodMT) / float(sumMT)
            else:
                mtCovDepthFinal = 0
            mtCovDepthFinal = "{0:.4f}".format(float(mtCovDepthFinal))
            if float(mtCovDepthFinal) < 0.0001:
                mtCovDepthFinal = "Below 0.0001"
            print("\tMT genome covered: " + str(mtCovered) + "%")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tMT genome covered: " + str(mtCovered) + "%")
            print("\tMT coverage depth: " + mtCovDepthFinal + "X")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\tMT coverage depth: " + mtCovDepthFinal + "X")
            coverages.write(
                sample_id + "\t" + genomeCovered + "%\t" + genCovDepthFinal + "X\t" + str(mtCovered) + "%\t" + mtCovDepthFinal + "X\n")


def admix(range_of_k_low, range_of_k_high, reorder=False):
    """
    Estimates ancestry in a model-based manner from large autossomal SNP datasets, using software ADMIXTURE.
    K-value must be specified. If one value only, cross-validation is automatically disabled.

    Input:
    Foo.bed
    (Foo.bim)
    (Foo.fam)


    Output:
    Foo.ancestry2.txt
    Foo.2.Q
    Foo.2.P

    :param range_of_k_low:
    (int)   Integer lower range value for K.

    :param range_of_k_high:
    (int)   Integer upper range value for K.

    :param reorder:
    (boolean)   FALSE if reordering of the individual to the end of the file is not wanted.
    (string)    Path to the text file containing the order of the individuals.
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< ADMIXTURE - Estimate ancestry >")
    fout.write("\n< ADMIXTURE - Estimate ancestry >")

    ## Checks if "bin4pca" folder already exists. Abortes if it doesn't
    if "bin4pca" in os.listdir(cwd):
        pass
    else:
        print("\t> ERROR: No 'bin4pca' folder found in directory, possible binary PLINK files inexistent. Aborting.")
        fout.write(
            "\n\t> ERROR: No 'bin4pca' folder found in directory, possible binary PLINK files inexistent. Aborting.")
        quit()

    ## Checks files and creates list of samples and list of files
    cwdBIN = cwd + "/bin4pca"
    samples_list = []
    files_list = []
    for i in os.listdir(cwdBIN):
        if "mergedex.bed" == i[-12:]:
            sample_id, sep, rest = i.partition("mergedex")
            if sample_id not in samples_list:
                samples_list.append(sample_id)
                files_list.append(i)
        elif "-merge.bed" == i[-10:]:
            sample_id, sep, rest = i.partition("-merge")
            if sample_id not in samples_list:
                samples_list.append(sample_id)
                files_list.append(i)

    ## If no input files, aborts here. Of they exist but folder "admix" does not, creates it
    if len(files_list) != 0 and "admix" not in os.listdir(cwdBIN):
        os.makedirs(cwdBIN + "/admix")
    elif len(samples_list) == 0:
        print("\t> ERROR: No input 'mergedex.bed' or '-merge.bed' files found in 'bin4pca' directory. Aborting.")
        fout.write("\n\t> ERROR: No input 'mergedex.bed' or '-merge.bed' files found in 'bin4pca' directory. Aborting.")
        quit()

    ## Reorder alphabetically and sepparate individuals in *.fam file to put ancient ones at the end.
    ## And places sample at the end of the list.
    if isinstance(reorder, str):
        for i in files_list:
            if "mergedex.bed" == i[-12:] or "-merge.bed" == i[-10:]:
                if "mergedex.bed" == i[-12:]:
                    sample_id, sep, rest = i.partition("mergedex")
                if "-merge.bed" == i[-10:]:
                    sample_id, sep, rest = i.partition("-merge")
                try:
                    test = open(cwdBIN + "/" + sample_id + sep + ".bim", "r")
                    test.close()
                    test2 = open(cwdBIN + "/" + sample_id + sep + ".fam", "r")
                    test2.close()
                    print("\t> Copying BED/BIM/FAM files into 'admix' folder...")
                    fout.write("\n\n\t> Copying BED/BIM/FAM files into 'admix' folder...")
                    shutil.copy2(cwdBIN + "/" + sample_id + sep + ".fam", cwdBIN + "/admix/" + sample_id + sep + ".fam")
                    shutil.copy2(cwdBIN + "/" + sample_id + sep + ".bim", cwdBIN + "/admix/" + sample_id + sep + ".bim")
                    shutil.copy2(cwdBIN + "/" + sample_id + sep + ".bed", cwdBIN + "/admix/" + sample_id + sep + ".bed")
                    shutil.copy2(reorder, cwdBIN + "/admix/" + sample_id + sep + "_reordered.fam")
                    fam_to_reorder = open(cwdBIN + "/" + sample_id + sep + ".fam", "r")
                    reordered_fam = open(cwdBIN + "/admix/" + sample_id + sep + "_reordered.fam", "a")
                    for line in fam_to_reorder:
                        if " " + sample_id + " " not in line:
                            pass
                        else:
                            fam, ind, a, b, c, d = line.split(" ")
                            our_line = fam + "\t" + ind + "\n"
                    reordered_fam.write(our_line)
                    reordered_fam.close()
                except FileNotFoundError:
                    print(
                        "\t> WARNING: '" + sample_id + sep + ".bim/fam' or 'param: reorder' file not found. Skipping " + sample_id + " admixture estimation.")
                    fout.write(
                        "\n\t> WARNING: '" + sample_id + sep + ".bim/fam' or 'param: reorder' file not found. Skipping " + sample_id + " admixture estimation.")

    ## Do not reorder the individuals in the *.fam file for situations where they come ordered and sepparated (e.g. anc from modern) from the dataset itself.
    ## Only places sample at the end of the list.
    elif reorder == False:
        for i in files_list:
            if "mergedex.bed" == i[-12:] or "-merge.bed" == i[-10:]:
                if "mergedex.bed" == i[-12:]:
                    sample_id, sep, rest = i.partition("mergedex")
                if "-merge.bed" == i[-10:]:
                    sample_id, sep, rest = i.partition("-merge")
                try:
                    test = open(cwdBIN + "/" + sample_id + sep + ".bim", "r")
                    test.close()
                    test2 = open(cwdBIN + "/" + sample_id + sep + ".fam", "r")
                    test2.close()
                    print("\t> Copying BED/BIM/FAM files into 'admix' folder...")
                    fout.write("\n\n\t> Copying BED/BIM/FAM files into 'admix' folder...")
                    shutil.copy2(cwdBIN + "/" + sample_id + sep + ".fam", cwdBIN + "/admix/" + sample_id + sep + ".fam")
                    shutil.copy2(cwdBIN + "/" + sample_id + sep + ".bim", cwdBIN + "/admix/" + sample_id + sep + ".bim")
                    shutil.copy2(cwdBIN + "/" + sample_id + sep + ".bed", cwdBIN + "/admix/" + sample_id + sep + ".bed")
                    fam_to_reorder = open(cwdBIN + "/" + sample_id + sep + ".fam", "r")
                    reordered_fam = open(cwdBIN + "/admix/" + sample_id + sep + "_reordered.fam", "w")
                    for line in fam_to_reorder:
                        if " " + sample_id + " " not in line:
                            fam, ind, a, b, c, d = line.split(" ")
                            reordered_fam.write(fam + "\t" + ind + "\n")
                        else:
                            fam, ind, a, b, c, d = line.split(" ")
                            our_line = fam + "\t" + ind + "\n"
                    reordered_fam.write(our_line)
                    reordered_fam.close()
                except FileNotFoundError:
                    print(
                        "\t> WARNING: '" + sample_id + sep + ".bim/fam' not found in current folder. Skipping " + sample_id + " admixture estimation.")
                    fout.write(
                        "\n\t> WARNING: '" + sample_id + sep + ".bim/fam' not found in current folder. Skipping " + sample_id + " admixture estimation.")

    admixfolder = cwdBIN + "/admix/"
    for i in os.listdir(cwdBIN + "/admix"):
        if "mergedex.bed" == i[-12:] or "-merge.bed" == i[-10:]:
            if "mergedex.bed" == i[-12:]:
                sample_id, sep, rest = i.partition("mergedex")
            if "-merge.bed" == i[-10:]:
                sample_id, sep, rest = i.partition("-merge")

            print("\t> Reordering or placing " + sample_id + " at the end of BED/BIM/FAM files...")
            fout.write("\n\n\t> Reordering or placing " + sample_id + " at the end of BED/BIM/FAM files...")

            test0 = sum(1 for line in open(admixfolder + sample_id + sep + "_reordered.fam"))
            test1 = sum(1 for line in open(admixfolder + sample_id + sep + ".fam"))
            if test0 != test1:
                print(
                    "\t\t> WARNING: Individuals' order file (*._reordered.fam) does not contain same number of individuals as in *.fam file. Skipping " + sample_id + " reordering.")
                fout.write(
                    "\n\t\t> WARNING: Individuals' order file (*._reordered.fam) does not contain same number of individuals as in *.fam file. Skipping " + sample_id + " reordering.")
                quit()
            else:
                print(
                    "plink --bfile " + cwdBIN + "/admix/" + sample_id + sep + " --indiv-sort f " + cwdBIN + "/admix/" + sample_id + sep + "_reordered.fam --make-bed --out " + cwdBIN + "/admix/" + sample_id + "_reordered --allow-no-sex")
                os.system(
                    "plink --bfile " + cwdBIN + "/admix/" + sample_id + sep + " --indiv-sort f " + cwdBIN + "/admix/" + sample_id + sep + "_reordered.fam --make-bed --out " + cwdBIN + "/admix/" + sample_id + "_reordered --allow-no-sex")
                fout.write(time.strftime(
                    "\n\t%H:%M:%S %Z ") + "\t" + "plink --bfile " + cwdBIN + "/admix/" + sample_id + sep + " --indiv-sort f " + cwdBIN + "/admix/" + sample_id + sep + "_reordered.txt --make-bed --out " + cwdBIN + "/admix/" + sample_id + "_reordered --allow-no-sex")

    ## Run ADMIXTURE software for all K-values in provided range
    for i in os.listdir(admixfolder):
        if "_reordered.bed" == i[-14:]:
            sample_id, sep, rest = i.partition("_reordered")
            print("\t> Running ADMIXTURE for " + sample_id + " with k-interval of " + str(range_of_k_low) + "-" + str(
                range_of_k_high))
            fout.write(
                "\n\n\t> Running ADMIXTURE for " + sample_id + " with k-interval of " + str(range_of_k_low) + "-" + str(
                    range_of_k_high))
            os.chdir(admixfolder)

            for s in range(range_of_k_low, range_of_k_high + 1):
                print("admixture --cv " + i + " " + str(s) + " > " + sample_id + "_cv" + str(s) + ".txt")
                os.system("admixture --cv " + i + " " + str(s) + " > " + sample_id + "_cv" + str(s) + ".txt")
                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "admixture --cv " + i + " " + str(
                    s) + " > " + sample_id + "_cv" + str(s) + ".txt")

    ## Open CV text files and scavenge K probability values
    print("\n\t> Collecting cross-validation error probabilities into 'Admix_CV.txt'")
    fout.write("\n\n\t> Collecting cross-validation error probabilities into 'Admix_CV.txt'")
    probs_list = {}
    samples_list2 = []
    for i in os.listdir(admixfolder):
        for s in range(range_of_k_low, range_of_k_high + 1):
            if "_cv" + str(s) + ".txt" in i:
                cv_file = open(i, "r")
                sample_id, sep, rest = i.partition("_cv")
                all_lines = cv_file.readlines()
                CV_line = [s for s in all_lines if "CV error (K=" in s]
                CV_line = CV_line[0].rstrip()
                sep, Kvalue, prob = re.split("CV error |: ", str(CV_line))
                Kvalue = Kvalue.strip("(")
                Kvalue = Kvalue.strip(")")
                probs_list[sample_id + ">>" + Kvalue] = prob
                if sample_id not in samples_list2:
                    samples_list2.append(sample_id)

    cv_out = open("Admix_CV" + ".txt", "w")
    for i in samples_list2:
        cv_out.write(
            "ADMIXTURE cross-validation errors for " + i + " for K interval of " + str(range_of_k_low) + "-" + str(
                range_of_k_high) + "\n")
        for s in range(range_of_k_low, range_of_k_high + 1):
            try:
                cv_out.write(i + ">>K=" + str(s) + "\t" + probs_list[i + ">>K=" + str(s)] + "\n")
            except KeyError:
                cv_out.write("No CV file for K=" + str(s) + " found for " + i + "." + "\n")
        cv_out.write("\n")

    cv_out.close()
    ## Can add in the future lines to automatically do the R plotting, but not for now


def admix2(range_of_k_low, range_of_k_high, multithread_cores, input_file="Laz16ALLPolarizeEIG1200KAutosomesOrdered"):
    """
    Estimates ancestry in a model-based manner from large autossomal SNP datasets, using software ADMIXTURE.
    K-value must be specified. If one value only, cross-validation is automatically disabled.

    Input:
    Foo.bed
    (Foo.bim)
    (Foo.fam)


    Output:
    Foo.ancestry2.txt
    Foo.2.Q
    Foo.2.P

    :param range_of_k_low:
    (int)   Integer lower range value for K.

    :param range_of_k_high:
    (int)   Integer upper range value for K.

    :param reorder:
    (boolean)   FALSE if reordering of the individual to the end of the file is not wanted.
    (string)    Path to the text file containing the order of the individuals.
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< ADMIXTURE - Estimate ancestry >")
    fout.write("\n< ADMIXTURE - Estimate ancestry >")

    ## Run ADMIXTURE software for all K-values in provided range
    print("\t> Running ADMIXTURE for " + input_file + " with k-interval of " + str(range_of_k_low) + "-" + str(
        range_of_k_high))
    fout.write("\n\n\t> Running ADMIXTURE for " + input_file + " with k-interval of " + str(range_of_k_low) + "-" + str(
        range_of_k_high))

    for s in range(range_of_k_low, range_of_k_high + 1):
        print("admixture --cv " + input_file + ".bed " + str(s) + " -j" + str(multithread_cores) + " > " + input_file + "_cv" + str(s) + ".txt")
        os.system("admixture --cv " + input_file + ".bed " + str(s) + " -j" + str(multithread_cores) + " > " + input_file + "_cv" + str(s) + ".txt")
        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "admixture --cv " + input_file + ".bed " + str(s) + " -j" + str(multithread_cores) + " > " + input_file + "_cv" + str(s) + ".txt")

    ## Open CV text files and scavenge K probability values
    print("\n\t> Collecting cross-validation error probabilities into 'Admix_CV.txt'")
    fout.write("\n\n\t> Collecting cross-validation error probabilities into 'Admix_CV.txt'")
    probs_list = {}
    samples_list2 = []
    for i in os.listdir(cwd):
        for s in range(range_of_k_low, range_of_k_high + 1):
            if "_cv" + str(s) + ".txt" in i:
                cv_file = open(i, "r")
                sample_id, sep, rest = i.partition("_cv")
                all_lines = cv_file.readlines()
                CV_line = [s for s in all_lines if "CV error (K=" in s]
                CV_line = CV_line[0].rstrip()
                sep, Kvalue, prob = re.split("CV error |: ", str(CV_line))
                Kvalue = Kvalue.strip("(")
                Kvalue = Kvalue.strip(")")
                probs_list[sample_id + ">>" + Kvalue] = prob
                if sample_id not in samples_list2:
                    samples_list2.append(sample_id)

    cv_out = open("Admix_CV" + ".txt", "w")
    for i in samples_list2:
        cv_out.write(
            "ADMIXTURE cross-validation errors for " + i + " for K interval of " + str(range_of_k_low) + "-" + str(
                range_of_k_high) + "\n")
        for s in range(range_of_k_low, range_of_k_high + 1):
            try:
                cv_out.write(i + ">>K=" + str(s) + "\t" + probs_list[i + ">>K=" + str(s)] + "\n")
            except KeyError:
                cv_out.write("No CV file for K=" + str(s) + " found for " + i + "." + "\n")
        cv_out.write("\n")

    cv_out.close()
    ## Can add in the future lines to automatically do the R plotting, but not for now


def admix3(k_num, reps, multithread_cores, input_file):
    """
    Estimates ancestry in a model-based manner from large autossomal SNP datasets, using software ADMIXTURE.

    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< ADMIXTURE - Estimate ancestry >")
    fout.write("\n< ADMIXTURE - Estimate ancestry >")

    ## Run ADMIXTURE software for all K-values in provided range
    print("\t> Running ", str(reps), " repetitions of ADMIXTURE for " + input_file + " with k of " + str(k_num))
    fout.write("\n\n\t> Running " + str(reps) + " repetitions of ADMIXTURE for " + input_file + " with k of " + str(k_num))

    for s in list(range(1, reps + 1)):
        print("admixture --cv " + input_file + ".bed " + str(k_num) + " -j" + str(multithread_cores) + " -s time > " + input_file + "_cv" + str(k_num) + "_" + str(s) + ".txt")
        os.system("admixture --cv " + input_file + ".bed " + str(k_num) + " -j" + str(multithread_cores) + " -s time > " + input_file + "_cv" + str(k_num) + "_" + str(s) + ".txt")
        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "admixture --cv " + input_file + ".bed " + str(k_num) + " -j" + str(multithread_cores) + " -s time > " + input_file + "_cv" + str(k_num) + "_" + str(s) + ".txt")
        os.rename(input_file + "." + str(k_num) + ".Q", input_file + "." + str(k_num) + ".r" + str(s) + ".Q")
        os.rename(input_file + "." + str(k_num) + ".P", input_file + "." + str(k_num) + ".r" + str(s) + ".P")


def end_analysis():
    """
    Outputs a tabbed file with sample name,number of total, trimmed and aligned reads, % of endogenous DNA, GC contents, mapDamage ratios and sex assessment.

    Input:
    *various_files*

    Output:
    Reads_Summary.txt

    """
    with open("Reads_Summary.txt", 'w') as foutr:
        os.chdir(originalWD)
        cwd = os.getcwd()
        list_samples = []
        print("\n< SUMMARY - Generating *.txt tabbed file with summary >")
        foutr.write(
            "Sample\tTotal_reads\t%GC\tTrimmed_reads\tAligned_reads\t%Endogenous\tUnique_Reads\t%Endogenous_Unique\t%GC\tMapDamage\tSex\t\n")
        # foutr.write("Sample\tTotal reads\tTrimmed reads\tAligned reads\t%End\t%GC\tAfter rmdup\t%Unique\t%GC\tMapDamage\tSex\t\n")
        fout.write("\n< SUMMARY - Generating *.txt tabbed file with summary >")
        warnings = []
        # for i in os.listdir(cwd):
        #     if 'fastq.gz' in i:
        #         sample_idd, sep, rest = i.partition('_L001')
        #         if sample_idd[-2:].isdigit() == True:
        #             sample_true_id = sample_idd[:-4]
        #         else:
        #             sample_true_id = sample_idd[:-3]
        #         list_samples.append(sample_true_id)

        for i in os.listdir(cwd):
            if '.trimmed.log' in i:
                sample_idd, sep, rest = i.partition('%')
            #if '.q30.bam' in i:
            #    sample_idd, sep, rest = i.partition('.q30.')
                list_samples.append(sample_idd)
                # list_samples.append(sample_true_id)

        for sample in list_samples:
            print("\n\t> Compiling data for sample " + sample)
            fout.write("\n\t> Compiling data for sample " + sample + "\n")
            samenametrimmedlog = 0
            warningslogfiles = []
            ## Processed and trimmed reads
            if any("trimmed.log" in i and sample in i for i in os.listdir(cwd)):
                for i in os.listdir(cwd):

                    if "trimmed.log" in i and sample in i:
                        sample_true, sep, rest = i.rpartition("%")
                        samenametrimmedlog += 1
                        if sample_true == sample:
                            if os.path.exists(cwd + "/" + i):
                                with open(i, 'r') as fin:
                                    ## Processed reads
                                    for line in fin:
                                        # print(line)
                                        if "Total reads processed:" in line:
                                            text = line.split()
                                            ba = text[3]
                                            ba = ba.replace(',', '')
                                            # a = fin.readline()
                                            # a = fin.readline()
                                            # a = fin.readline()
                                            # a = fin.readline()
                                            # a = fin.readline()
                                            # ba = a[21:].strip()
                                            # print(ba)
                                            # print("Processed reads: "+ba)
                                            ## Trimmed reads
                                        if "Reads with adapters:" in line:
                                            # a = fin.readline()
                                            # a = fin.readline()
                                            # ca = line[21:32].strip()
                                            # print(line)
                                            text = line.split()
                                            ca = text[3]
                                            ca = ca.replace(',', '')
                                            # print("Trimmed reads: "+ca)
                        else:
                            warningslogfiles.append(
                                sample + " - No *.trimmed.log file for " + sample + " > No 'Processed reads' & 'Trimmed reads'")
                            ba = "N/A"
                            ca = "N/A"
                            pass

                if samenametrimmedlog <= len(warningslogfiles):
                    print(
                        "\t> WARNING: No *.trimmed.log file for " + sample + " in '" + cwd + "'. 'Processed reads' & 'Trimmed reads' can not be retrieved.")
                    warnings.append(
                        sample + " - No *.trimmed.log file for " + sample + " > No 'Processed reads' & 'Trimmed reads'")
                else:
                    pass
            else:
                print(
                    "\t> WARNING: No *.trimmed.log file for " + sample + " in '" + cwd + "'. 'Processed reads' & 'Trimmed reads' can not be retrieved.")
                warnings.append(
                    sample + " - No *.trimmed.log file for " + sample + " > No 'Processed reads' & 'Trimmed reads'")
                ba = "N/A"
                ca = "N/A"

            ## Aligned reads
            if any("sort_fastqc.html" in i and sample in i for i in os.listdir(cwd)):
                for i in os.listdir(cwd):
                    if "sort_fastqc.html" in i and sample in i:
                        with open(i, 'r') as fin6:
                            ## Aligned reads
                            for line in fin6:
                                if "Total Sequences</td><td>" in line:
                                    # text = line.split()
                                    # print(line.find("Total Sequences</td><td>"))
                                    text = (line[line.find("Total Sequences</td><td>"):line.find(
                                        "Total Sequences</td><td>") + 26 + 60])
                                    # print(text)
                                    # print(text.split("</td><td>"))
                                    reads = text.split("</td><td>")[1].split("</td></tr><tr><td>")[0]
                                    reads = reads.replace(',', '')
                                    # ba = text[3]
            else:
                print("\t> WARNING: No FASTQC *.sort file in '" + cwd + "'. 'Aligned reads' can not be calculated.")
                warnings.append(sample + " - No FASTQC *.sort file for " + sample + " > No 'aligned reads'")
                reads = "N/A"
                pass

            ## Reads after duplicates removal and GC content
            if any("rmdup_fastqc.html" in i and sample in i for i in os.listdir(cwd)):
                for i in os.listdir(cwd):
                    if "rmdup_fastqc.html" in i and sample in i:
                        with open(i, 'r') as fin3:
                            ## Reads after rmdup
                            for line in fin3:
                                if "Total Sequences</td><td>" in line:
                                    text = (line[line.find("Total Sequences</td><td>"):line.find(
                                        "Total Sequences</td><td>") + 26 + 60])
                                    # print(text)
                                    # print(text.split("</td><td>"))
                                    rmdpreads = text.split("</td><td>")[1].split("</td></tr><tr><td>")[0]
                                    rmdpreads = rmdpreads.replace(',', '')
                                    text2 = (line[line.find("%GC</td><td>"):line.find("%GC</td><td>") + 13 + 60])
                                    gc_cont_rmdup = text2.split("</td><td>")[1].split("</td></tr>")[0]
            else:
                print("\t> WARNING: No FASTQC *.rmdup file in '" + cwd + "'. 'Unique reads' can not be calculated.")
                warnings.append(sample + " - No FASTQC *.rmdup file for " + sample + " > No 'unique reads'")
                rmdpreads = "N/A"
                gc_cont_rmdup = "N/A"
                pass

            # ## Initial GC content
            if any("_fastqc.html" in i and sample in i and "sort" not in i and "rmdup" not in i for i in
                   os.listdir(cwd)):
                for i in os.listdir(cwd):
                    if "_fastqc.html" in i and sample in i and "sort" not in i and "rmdup" not in i:
                        with open(i, 'r') as fin4:
                            for line in fin4:
                                # print(line)
                                if "Total Sequences</td><td>" in line:
                                    # print(text)
                                    # print(text.split("</td><td>"))
                                    text2 = (line[line.find("%GC</td><td>"):line.find("%GC</td><td>") + 13 + 60])
                                    gc_cont = text2.split("</td><td>")[1].split("</td></tr>")[0]
                                    # print("Initial GC cont: "+gc_cont)
            else:
                print(
                    "\t> WARNING: No FASTQC initial file in '" + cwd + "'. 'Initial GC content' can not be retrieved.")
                warnings.append(sample + " - No FASTQC initial file for " + sample + " > No 'Initial GC content'")
                gc_cont = "N/A"
                pass

            ## Endogenous contents by total reads
            if reads != "N/A" and ba != "N/A" and rmdpreads != "N/A":
                # if reads != "N/A" and ca != "N/A" and rmdpreads != "N/A":
                perc_end = (int(reads) / int(ba)) * 100
                # print("% End:","%.2f" % perc_end)

                perc_uni = (int(rmdpreads) / int(ba)) * 100
                # print("% Unique:","%.2f" % perc_uni)
            else:
                perc_end = "N/A"
                perc_uni = "N/A"
                # print("% End: N/A\n% Unique: N/A")

            ## Mapdamage 5' C>T and 3' G>A ratios
            if os.path.exists(cwd + "/mapdamage/mapdamage_" + sample):
                file10 = cwd + "/mapdamage/mapdamage_" + sample + "/5pCtoT_freq.txt"
                with open(file10, 'r') as fin10:
                    positions = [1]
                    fin10.readline()
                    for line in fin10:
                        position, sep, freq = line.partition("\t")
                        freq = freq.strip("\n")
                        positions.append(freq)
                    if positions[1] > positions[2] > positions[3] and positions[1] > positions[25] and all(
                                    x != "0" for x in positions):
                        positionssm = positions[4:len(positions)]
                        if all(i < positions[3] for i in positionssm) == True:
                            #map5ctratio = float(positions[1]) - float(positions[25])
                            map5ctratio = float(positions[1])
                            print(
                                "\t> Ancient DNA damage pattern potentially detected. 5' C>T transition differential bp1-bp25: " + "%.2f" % float(
                                    map5ctratio))
                            fout.write(time.strftime(
                                "\t%H:%M:%S %Z ") + "\tAncient DNA damage pattern potentially detected. 5' C>T transition differential bp1-bp25: " + "%.2f" % float(
                                map5ctratio))
                        else:
                            print("\t> Ancient DNA damage pattern NOT found for 5' C>T. Please confirm.")
                            map5ctratio = "N/A"
                    else:
                        print("\t> Ancient DNA damage pattern NOT found for 5' C>T. Please confirm.")
                        map5ctratio = "N/A"

                file11 = cwd + "/mapdamage/mapdamage_" + sample + "/3pGtoA_freq.txt"
                with open(file11, 'r') as fin11:
                    positions3 = [1]
                    fin11.readline()
                    for line in fin11:
                        position, sep, freq = line.partition("\t")
                        freq = freq.strip("\n")
                        positions3.append(freq)
                    if positions3[1] > positions3[2] > positions3[3] and positions3[1] > positions3[25] and all(
                                    z != "0" for z in positions3):
                        positionssm3 = positions3[4:len(positions3)]
                        if all(i < positions[3] for i in positionssm3) == True:
                            #map3garatio = float(positions3[1]) - float(positions3[25])
                            map3garatio = float(positions3[1])
                            print(
                                "\t> Ancient DNA damage pattern potentially detected. 3' G>A transition differential bp1-bp25: " + "%.2f" % float(
                                    map3garatio))
                            fout.write(time.strftime(
                                "\n\t%H:%M:%S %Z ") + "\tAncient DNA damage pattern potentially detected. 3' G>A transition differential bp1-bp25: " + "%.2f" % float(
                                map3garatio) + "\n")
                        else:
                            print("\t> Ancient DNA damage pattern NOT found for 3' A>G. Please confirm.")
                            map3garatio = "N/A"
                    else:
                        print("\t> Ancient DNA damage pattern NOT found for 3' A>G. Please confirm.")
                        map3garatio = "N/A"
            else:
                map5ctratio = "N/A"
                map3garatio = "N/A"
                print("\t> WARNING: No mapDamage folder detected for " + sample + ". Please confirm.")
                warnings.append(sample + " - No mapDamage folder detected for " + sample + ". Please confirm.")
                pass

            foutr.write(sample)
            foutr.write("\t")
            foutr.write(ba)
            foutr.write("\t")
            foutr.write(gc_cont)
            foutr.write("\t")
            foutr.write(ca)
            foutr.write("\t")
            foutr.write(reads)
            foutr.write("\t")
            if perc_end != "N/A":
                foutr.write("%.2f" % perc_end)
            else:
                foutr.write("N/A")
            foutr.write("\t")
            foutr.write(rmdpreads)
            foutr.write("\t")
            if perc_uni != "N/A":
                foutr.write("%.2f" % perc_uni)
            else:
                foutr.write("N/A")
            foutr.write("\t")
            foutr.write(gc_cont_rmdup)
            foutr.write("\t")
            if map5ctratio != "N/A" and map3garatio != "N/A":
                foutr.write("%.2f" % float(map5ctratio))
            else:
                foutr.write("N/A")
            foutr.write(" | ")
            if map3garatio != "N/A" and map5ctratio != "N/A":
                foutr.write("%.2f" % float(map3garatio))
            else:
                foutr.write("N/A")
            foutr.write("\t")

            ## Sexing (Skoglund 2014)
            if os.path.exists(cwd + "/sexing_" + sample + ".txt"):
                sexfile = cwd + "/sexing_" + sample + ".txt"

                with open(sexfile, 'r') as fin5:
                    try:
                        if os.stat(sexfile).st_size == 0:
                            # print("Sexing: File '"+i+"' is empty")
                            foutr.write("-")
                        else:
                            s = fin5.readlines()[1]
                            slist = s.split("\t")
                            s = slist[6].strip("\n")
                            ychr = slist[2].strip(" ")
                            yxchr = slist[1].strip(" ")
                            ci = slist[5].strip(" ")
                            if ci != "0.0-0.0":
                                # print("Sexing: "+s[0].capitalize()+s[1:]+", with 95% CI:"+ci+" (chrY:"+ychr+"/"+yxchr+")")
                                foutr.write(s[0].capitalize())
                                foutr.write(s[1:])
                                foutr.write(", with 95% CI:")
                                foutr.write(ci)
                                foutr.write(" (chrY:")
                                foutr.write(ychr)
                                foutr.write("/")
                                foutr.write(yxchr)
                                foutr.write(")")
                                # foutr.write("\n")
                            else:
                                # print("Sexing: Inconclusive, with 95% CI:"+ci+", but "+s+" (chrY:"+ychr+"/"+yxchr+")")
                                foutr.write("Inconclusive")
                                foutr.write(", with 95% CI:")
                                foutr.write(ci)
                                foutr.write(", but ")
                                foutr.write(s)
                                foutr.write(" (chrY:")
                                foutr.write(ychr)
                                foutr.write("/")
                                foutr.write(yxchr)
                                foutr.write(")")
                                # foutr.write("\n")
                    except IndexError:
                        print("\t> WARNING: Incorrect sexing file for " + sample + ". No sexual assignment.")
                        warnings.append(sample + " - Incorrect sexing file for " + sample + ". No sexual assignment.")
                        foutr.write("N/A")
                        pass

            elif os.path.exists(cwd + "/Sexing_stats.txt"):
                sexfile = cwd + "/Sexing_stats.txt"
                with open(sexfile, 'r') as fin5:
                    try:
                        if os.stat(sexfile).st_size == 0:
                            # print("Sexing: File '"+i+"' is empty")
                            foutr.write("-")
                        else:
                            for line in fin5:
                                if "Based on Jones" not in line:
                                    # print(line)
                                    samp, xx, xy, pv, assessxy, assess = line.split("\t")
                                    if samp == sample:
                                        if float(pv) <= 0.05:
                                            foutr.write(assess.rstrip("\n") + " with p-value of " + pv)
                    except IndexError:
                        print("\t> WARNING: Incorrect sexing file for " + sample + ". No sexual assignment.")
                        warnings.append(sample + " - Incorrect sexing file for " + sample + ". No sexual assignment.")
                        foutr.write("N/A")
                        pass

            else:
                print("\t> WARNING: No sexing file for " + sample + ". No sexual assignment.")
                warnings.append(sample + " - No sexing file for " + sample + ". No sexual assignment.")
                foutr.write("N/A")
                pass

            foutr.write("\n")

        if len(warnings) > 0:
            print("\n\t< WARNINGS - Please review your data >")
            warnings = sorted(list(set(warnings)))
            fout.write("\n\t< WARNINGS Please review your data >")
            for i in warnings:
                print("\t\t>", i)

                fout.write("\n\t\t> " + i)


def subsampling(fract, iterations, forSexing=False):  #### MAKE IT PRETTIER

    cwd = os.getcwd()

    for i in os.listdir(cwd):
        if "rmdup.bam" in i:
            #sampID, sep, rest = i.rpartition(".q30")
            sampID, rest = re.split('\.q..', i)
            sep = re.findall('\.q..', i)[0]
            zero, sep2, fractsep = str(fract).rpartition(".")
            seeds = list(range(500000))
            #print(sep)
            #print(seeds)
            #print(i)
            #print(sampID)
            seedsrand = (random.sample(seeds, iterations))
            #print(seedsrand)
            it = 1
            for ite in seedsrand:
                print("samtools view -s " + str(ite) + "." + str(fractsep) + " -b " + i + " > " + sampID + sep + ".rmdup.sub" + str(it) + ".bam")
                os.system("samtools view -s " + str(ite) + "." + str(fractsep) + " -b " + i + " > " + sampID + sep + ".rmdup.sub" + str(it) + ".bam")
                it += 1

                if forSexing == True:
                    for file in os.listdir(cwd):
                        if sampID + sep + ".rmdup.sub" + str(it - 1) + ".bam" in file:
                            print('\tsamtools index ' + file + ' ' + file + '.bai')
                            os.system(r'samtools index ' + file + ' ' + file + '.bai')
                            print("\t\t" + r"samtools idxstats " + file + " > " + sampID + ".sub" + str(it - 1) + ".chr.txt\n")
                            os.system(r"samtools idxstats " + file + " > " + sampID + ".sub" + str(it - 1) + ".chr.txt")
                            #                    fout.write("\n\t> Reading distribution of reads per chromosome from '"+i+"'..."+time.strftime("\n\t%H:%M:%S %Z ")+"\t"+r"samtools idxstats "+i+" > "+sample_id+".chr.txt\n")

                elif forSexing == False:
                    pass


def extract_reads_pos(chrom, positionL, positionH, outputBam=False):
    """
    Extracts reads that contain the specific genomic position or range provided. Repeat lower and higher positions input to extract a single position.

    Input:
    Foo.bam

    Output:
    Foo.2_123455-123455.bam

    :param chrom:
    (int)
    (str)   Integer or string (for mitochondrial DNA) of the chromosome number.

    :param positionL:
    (int)   Integer of the desired contig lower 5' position.

    :param positionH:
    (int)   Integer of the desired contig higher 3' position.
    """

    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< EXTRACT READS POS - Extracting reads aligning to specific genome positions >")
    fout.write("\n< EXTRACT READS POS - Extracting reads aligning to specific genome positions >")
    print("\t> Extracting reads containing position " + str(chrom) + ":" + str(positionL) + "-" + str(positionH))
    fout.write("\t> Extracting reads containing position " + str(chrom) + ":" + str(positionL) + "-" + str(positionH))

    for i in os.listdir(cwd):
        if "_merged" in i:
            if re.search("_merged\.q..\.rmdup\.bam$", i) != None:
            #if "_merged.q30.rmdup.bam" == i[-21:]:
                sample_id, sep, rest = i.partition("_merged.")
                sep = re.findall('\.q..', i)[0]
                print("\t> Extracting reads for " + sample_id)
                fout.write("\n\n\t> Extracting reads for " + sample_id)
                print("\tsamtools view -b " + i + " " + "chr" + str(chrom) + ":" + str(positionL) + "-" + str(positionH) + " > " + sample_id + ".filter." + str(chrom) + "_" + str(positionL) + "-" + str(positionH) + sep + ".rmdup.bam ")
                os.system("samtools view -b " + i + " " + "chr" + str(chrom) + ":" + str(positionL) + "-" + str(positionH) + " > " + sample_id + ".filter." + str(chrom) + "_" + str(positionL) + "-" + str(positionH) + sep + ".rmdup.bam ")
                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools view -b " + i + " " + "chr" + str(chrom) + ":" + str(positionL) + "-" + str(positionH) + " > " + sample_id + ".filter." + str(chrom) + "_" + str(positionL) + "-" + str(positionH) + sep + ".rmdup.bam ")
            else:
                pass
        elif re.search("\.q..\.rmdup\.bam$", i) != None and "merge" not in i:
            #sample_id, sep, rest = i.partition(".q30")
            sample_id, rest = re.split('\.q..', i)
            sep = re.findall('\.q..', i)[0]
            print("\t> Extracting reads for " + sample_id)
            fout.write("\n\n\t> Extracting reads for " + sample_id)
            print("\tsamtools view -b " + i + " " + "chr" + str(chrom) + ":" + str(positionL) + "-" + str(positionH) + " > " + sample_id + ".filter." + str(chrom) + "_" + str(positionL) + "-" + str(positionH) + sep + ".rmdup.bam ")
            os.system("samtools view -b " + i + " " + "chr" + str(chrom) + ":" + str(positionL) + "-" + str(positionH) + " > " + sample_id + ".filter." + str(chrom) + "_" + str(positionL) + "-" + str(positionH) + sep + ".rmdup.bam ")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools view -b " + i + " " + "chr" + str(chrom) + ":" + str(positionL) + "-" + str(positionH) + " > " + sample_id + ".filter." + str(chrom) + "_" + str(positionL) + "-" + str(positionH) + sep + ".rmdup.bam ")


def extract_reads_seq(sequ, thres, mismatchBP, oneFile=False):
    """
    Extracts reads that contain the specific nucleotide sequence provided.
    NOTE: Needs "regex" module (pypi.python.org/pypi/regex/).

    Input:
    Foo.fastq

    Output:
    aaaaaaaaaaaaa

    :param sequ:
    (str)   String with DNA sequence to search for.

    :param oneFile:
    (boolean)
    (str)   String with name of single file to check. Checks all files in folder if "False".

    :param thres:
    (int)   Integer with minimum match threshold (as fraction).
    """

    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< EXTRACT READS SEQ - Extracting reads containing specific sequence >")
    fout.write("\n< EXTRACT READS SEQ - Extracting reads containing specific sequence >")
    print("\t> Extracting reads containing sequence " + str(sequ))
    fout.write("\n\t> Extracting reads containing sequence " + str(sequ))
    print("\t> Minimum match threshold: " + str(thres))
    fout.write("\n\t> Minimum match threshold: " + str(thres))
    print("\t> Allowed base-pair mismatches: " + str(mismatchBP))
    fout.write("\n\t> Allowed base-pair mismatches: " + str(mismatchBP))
    counthits = 0
    countPartialHits = 0
    counter = 1
    lineCounter = 1

    ## Defines a function that returns all the possible sequences of letters in a given string
    def substr(string):
        j = 1
        a = set()
        while True:
            for i in range(len(string) - j + 1):
                a.add(string[i:i + j])
            if j == len(string):
                break
            j += 1
        return a

    ## Creates lists of possible substrings of the provided DNA sequence for the 5' and 3' end based on the mininum match threshold
    ## For sequence 'ATCTAGCGGATG' it will create the following lists with a 0.5 threshold:
    ## 5' - 'ATCTAGCGGATG','ATCTAGCGGAT','ATCTAGCGGA','ATCTAGCGG','ATCTAGCG','ATCTAGC','ATCTAG'
    ## 3' - 'ATCTAGCGGATG','TCTAGCGGATG','CTAGCGGATG','TAGCGGATG','AGCGGATG','GCGGATG','CGGATG'
    permuts = list(substr(sequ))
    permutsThres5 = list()
    permutsThres3 = list()
    minLen = int(len(sequ) * thres)
    for i in permuts:
        if sequ[0:minLen] in i and len(i) >= minLen:
            permutsThres5.append(i)
    for i in permuts:
        if sequ[minLen:] in i and len(i) >= minLen:
            permutsThres3.append(i)

    ## Goes through each read of given file
    if oneFile != False:
        print("\t> Searching sequence in file " + oneFile)
        fout.write("\n\t> Searching sequence in file " + oneFile)
        fastqin = open(oneFile, 'r')
        fastqlines = fastqin.readlines()
        fastqout = list()
        # print(len(fastqlines))
        for line in fastqlines:
            if line[0] == "@":
                # print("----------------------------")
                # print("Counter: ",counter)
                # print("Line: ", line)
                # print("FastqLine[5]: "+fastqlines[5])
                # print("Fastqline:", fastqlines[counter-1])
                # print("Fastqline+1:", fastqlines[counter])
                # print("Fastqline+2:", fastqlines[counter+1])
                # counter +=4
                #                print(counter)
                #                print(round(len(fastqlines)*0.01-4))
                #                 print(lineCounter)
                #                 print(counter)
                #                 print(round((len(fastqlines)/4)*0.2))
                if lineCounter == round((len(fastqlines) / 4) * 0.2):
                    print("\t\tLines processed: " + str(counter) + " of " + str(len(fastqlines)) + " (20%)")
                if lineCounter == round((len(fastqlines) / 4) * 0.4):
                    print("\t\tLines processed: " + str(counter) + " of " + str(len(fastqlines)) + " (40%)")
                if lineCounter == round((len(fastqlines) / 4) * 0.6):
                    print("\t\tLines processed: " + str(counter) + " of " + str(len(fastqlines)) + " (60%)")
                if lineCounter == round((len(fastqlines) / 4) * 0.8):
                    print("\t\tLines processed: " + str(counter) + " of " + str(len(fastqlines)) + " (80%)")
                lineCounter += 1
                inLineCounter = 0
                inLinePartialCounter = 0
                if counter <= len(fastqlines):
                    ja = regex.findall("(" + sequ + "){e<=" + str(mismatchBP) + "}", fastqlines[counter])
                    if len(ja) != 0:
                        # print(ja[0])
                        inLineCounter += 1
                        # input("Found COMPLETE")

                    if inLineCounter == 0:
                        for i5 in permutsThres5:
                            lengto = int(len(i5)) + 1
                            # print("i5: "+i5)
                            # print("Normal Line: "+fastqlines[counter])
                            # print("Short Line: "+fastqlines[counter][-lengto:])
                            ja5 = regex.findall("(" + i5 + "){e<=" + str(mismatchBP) + "}",
                                                fastqlines[counter][-lengto:])  # allows up to X base mismatches
                            # if i5 == fastqlines[counter+1][-lengto:]:
                            if len(ja5) != 0:
                                # print(ja5[0])
                                inLinePartialCounter += 1
                                # input("Found I5")

                        for i3 in permutsThres3:
                            lengto = int(len(i3))
                            # print("i3: "+i3)
                            # print("Normal line:"+fastqlines[counter])
                            # print("Short line:"+fastqlines[counter][0:lengto])
                            ja3 = regex.findall("(" + i3 + "){e<=" + str(mismatchBP) + "}",
                                                fastqlines[counter][0:lengto])  # allows up to X base mismatches
                            # if i3 == fastqlines[counter+1][0:lengto]:
                            if len(ja3) != 0:
                                # print(ja3[0])
                                # input("Found I3")
                                # counthits +=1
                                # countPartialHits +=1
                                inLinePartialCounter += 1

                if inLineCounter != 0:
                    fastqout.append(fastqlines[counter - 1])
                    fastqout.append(fastqlines[counter])
                    fastqout.append(fastqlines[counter + 1])
                    fastqout.append(fastqlines[counter + 2])
                    counthits += 1
                if inLineCounter == 0 and inLinePartialCounter != 0:
                    fastqout.append(fastqlines[counter - 1])
                    fastqout.append(fastqlines[counter])
                    fastqout.append(fastqlines[counter + 1])
                    fastqout.append(fastqlines[counter + 2])
                    countPartialHits += 1
                counter += 4

    print("\t\tNumber of hits: " + str(counthits))
    print("\t\tNumber of partial hits: " + str(countPartialHits))
    fout.write("\n\t\tNumber of hits: " + str(counthits))
    fout.write("\n\t\tNumber of partial hits: " + str(countPartialHits))
    if len(fastqout) != 0:
        fastqoutfile = open(oneFile + "_seq_search", 'w')
        for i in fastqout:
            fastqoutfile.write(i)


def hirisplex_lct(index, list_file="embedded"):
    """
    Extracts hirisplex system and lactose intolerance genotypes.

    Input:
    Foo.q30.rmdup.bam

    Output:
    Foo.pileupgatk.hirisplexlct.txt

    :param index:
    (str)   Path to reference genome fasta file.

    :param list_file:
    (str)   Deafult: "embedded", which uses an embedded list of positions. You can also provide the path to a text file containing list of the chromosomal locations of the genotypes to extract.
    """

    cwd = os.getcwd()
    embedded_list = ["chr16:89985750-89985754", "chr16:89986091-89986091", "chr16:89986154-89986154",
                     "chr16:89986144-89986144", "chr16:89985844-89985844", "chr16:89985918-89985918",
                     "chr16:89986117-89986117", "chr16:89986546-89986546", "chr16:89986122-89986122",
                     "chr16:89985940-89985940", "chr16:89986130-89986130", "chr5:33958959-33958959",
                     "chr5:33951693-33951693", "chr12:89328335-89328335", "chr6:457748-457748", "chr6:396321-396321",
                     "chr11:88911696-88911696", "chr15:28230318-28230318", "chr14:92801203-92801203",
                     "chr15:28365618-28365618", "chr20:33218090-33218090", "chr14:92773663-92773663",
                     "chr11:89011046-89011046", "chr9:12709305-12709305", "chr2:136608643-136608643",
                     "chr2:136608646-136608646", "chr2:136608651-136608651", "chr2:136616754-136616754",
                     "chr2:136608746-136608746", "chr20:32785212-32785212", "chr15:28187772-28187772",
                     "chr15:48426484-48426484"]
    if list_file == "embedded":
        hiris = open("HIRISPLEX_LCT.list", "w")
        for line in embedded_list:
            hiris.write(line + "\n")
        hiris.close()
        list_file = cwd + "/HIRISPLEX_LCT.list"

    for i in os.listdir(cwd):

        if re.search("\.q..\.rmdup\.bam$", i) != None:
        #if '.q30.rmdup.bam' == i[-14:]:
            #sample_id, sep, rest = i.partition('.q30')
            sample_id, rest = re.split('\.q..', i)
            sep = re.findall('\.q..', i)[0]
            print("\n\t> Collecting hirisplex PILEUP data from sample: " + sample_id)
            if i + '.bai' in os.listdir(cwd):
                print(
                    '\t\tgatk \n\tPileup \n\t-R ' + index + ' \n\t-L ' + list_file + ' \n\t-I ' + i + ' \n\t-O ' + sample_id + '.pileupgatk.hirisplexlct.txt')
                os.system(
                    'gatk Pileup -R ' + index + ' -L ' + list_file + ' -I ' + i + ' -O ' + sample_id + '.pileupgatk.hirisplexlct.txt')
                fout.write(time.strftime(
                    "\n\t%H:%M:%S %Z ") + 'gatk Pileup -R ' + index + ' -L ' + list_file + ' -I ' + i + ' -O ' + sample_id + '.pileupgatk.hirisplexlct.txt')


def hirisplexs(referenceGenome, minBQ = 30, minMQ = 30, exclTerminals = False, hpsconvertonline_location = "/home/dfernandes/Software/hpsconvertsonline.R"):
    """
    Extracts HIrisPlex-S genotypes and creates a CSV for direct upload to https://hirisplex.erasmusmc.nl/.
    Requires "hpsconvertonline.R" from https://walshlab.sitehost.iu.edu/pages/tools.html, and a few edits to it:
    
    Replace line 5:
    SourceFile<-file.choose()
    
    By:
    args = commandArgs(trailingOnly=TRUE)
    SourceFile = args[1]
    options(error=expression(NULL))

    And add these lines right before the last line with write.table:
    online[online == "NNANNA"] = 'NA'
    online[online == "NANA"] = 'NA'

    Input:
    Foo.q30.rmdup.bam

    Output:
    Foo.pileuphirisplexs.txt

    :param index:
    (str)   Path to reference genome fasta file.

    :param list_file:
    (str)   Deafult: "embedded", which uses an embedded list of positions. You can also provide the path to a text file containing list of the chromosomal locations of the genotypes to extract in BED format.
    """

    cwd = os.getcwd()
    import csv
    #embedded_list = ["5 33951692 33951693","5 33958958 33958959","6 396320 396321","6 457747 457748","9 12709304 12709305","9 16858083 16858084","11 88911695 88911696","11 89011045 89011046","11 89017960 89017961","12 89328334 89328335","14 92773662 92773663","14 92801202 92801203","14 92882825 92882826","15 28187771 28187772","15 28197036 28197037","15 28230317 28230318","15 28271774 28271775","15 28288120 28288121","15 28356858 28356859","15 28365617 28365618","15 28453214 28453215","15 28496194 28496195","15 28530181 28530182","15 48426483 48426484","16 89383724 89383725","16 89984377 89984378","16 89985843 89985844","16 89985917 89985918","16 89985939 89985940","16 89986090 89986091","16 89986116 89986117","16 89986121 89986122","16 89986129 89986130","16 89986143 89986144","16 89986153 89986154","16 89986545 89986546","16 90024205 90024206","2 32665747 32665748","2 32785211 32785212","2 33218089 33218090"]
    embedded_list = ['2 32665747 32665748','2 32785211 32785212','2 33218089 33218090','5 33951692 33951693','5 33958958 33958959','6 396320 396321','6 457747 457748','9 12709304 12709305','9 16858083 16858084','11 88911695 88911696','11 89011045 89011046','11 89017960 89017961','12 89328334 89328335','14 92773662 92773663','14 92801202 92801203','14 92882825 92882826','15 28187771 28187772','15 28197036 28197037','15 28230317 28230318','15 28271774 28271775','15 28288120 28288121','15 28356858 28356859','15 28365617 28365618','15 28453214 28453215','15 28496194 28496195','15 28530181 28530182','15 48426483 48426484','16 89383724 89383725','16 89984377 89984378','16 89985843 89985844','16 89985917 89985918','16 89985939 89985940','16 89986090 89986091','16 89986116 89986117','16 89986121 89986122','16 89986129 89986130','16 89986143 89986144','16 89986153 89986154','16 89986545 89986546','16 90024205 90024206']
    #embedded_list_rs_for_upload = ['rs16891982_C','rs28777_C','rs12203592_T','rs4959270_A','rs683_G','rs10756819_G','rs1042602_T','rs1393350_T','rs1126809_A','rs12821256_G','rs12896399_T','rs2402130_G','rs17128291_C','rs1545397_T','rs1800414_C','rs1800407_A','rs12441727_A','rs1470608_A','rs1129038_G','rs12913832_T','rs2238289_C','rs6497292_C','rs1667394_C','rs1426654_G','rs3114908_T','rs3212355_A','rs1805005_T','rs1805006_A','rs2228479_A','rs11547464_A','rs1805007_T','rs201326893_A','rs1110400_C','rs1805008_T','rs885479_T','rs1805009_C','rs8051733_C','rs6059655_T','rs6119471_C','rs2378249_C']
    embedded_list_rs_for_upload = ['rs6059655_T','rs6119471_C','rs2378249_C','rs16891982_C','rs28777_C','rs12203592_T','rs4959270_A','rs683_G','rs10756819_G','rs1042602_T','rs1393350_T','rs1126809_A','rs12821256_G','rs12896399_T','rs2402130_G','rs17128291_C','rs1545397_T','rs1800414_C','rs1800407_A','rs12441727_A','rs1470608_A','rs1129038_G','rs12913832_T','rs2238289_C','rs6497292_C','rs1667394_C','rs1426654_G','rs3114908_T','rs3212355_A','rs1805005_T','rs1805006_A','rs2228479_A','rs11547464_A','rs1805007_T','rs201326893_A','rs1110400_C','rs1805008_T','rs885479_T','rs1805009_C','rs8051733_C']

    hiris = open("hirisplexs.bed", "w")
    for line in embedded_list:
        hiris.write(line + "\n")
    hiris.close()
    list_file = cwd + "/hirisplexs.bed"

    print("\t< HIRISPLEX-S - Extract HIrisPlex-S genotypes and prepare file for upload >")
    fout.write("\n\n\t< HIRISPLEX-S - Extract HIrisPlex-S genotypes and prepare file for upload >")

    for i in os.listdir(cwd):
        if re.search("\.q..\.rmdup\.bam$", i) != None:
            sample_id, rest = re.split('\.q..', i)
            sep = re.findall('\.q..', i)[0]
            print("\n\t> Collecting HIrisPlex-S PILEUP data from sample: " + sample_id)
            if i + '.bai' in os.listdir(cwd):
                print("\tsamtools mpileup -a -Q " + str(minBQ) + " -q " + str(minMQ) + " -B -f " + referenceGenome + " -l " + list_file + " " + i + " > " + sample_id + ".pileuphirisplexs.txt")
                #os.system("samtools mpileup -a -Q " + str(minBQ) + " -q " + str(minMQ) + " -B -f " + referenceGenome + " -l " + list_file + " " + i + " > " + sample_id + ".pileuphirisplexs.txt")
                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools mpileup -a -Q " + str(minBQ) + " -q " + str(minMQ) + " -B -f " + referenceGenome + " -l " + list_file + " " + i + " > " + sample_id + ".pileuphirisplexs.txt")

    csvLine1 =["sampleid"]
    csvLine1.extend(embedded_list_rs_for_upload)
    #print(csvComb)
    #print(jarrrrr)
    csvOut = open('HIrisPlex-S_toConvert.csv', 'w+', newline = '')
    writer = csv.writer(csvOut)
    writer.writerow(csvLine1)
    
    for i in os.listdir(cwd):
        if '.pileuphirisplexs.txt' == i[-21:]:
            if os.path.getsize(i) > 0:
                sample_id, rest = re.split('.pileuphirisplexs.txt', i)

#                print("awk '{print $1" + '"-"' + "$2}' " + i + " > temporary1")
#                os.system("awk '{print $1" + '"-"' + "$2}' " + i + " > temporary1")
#                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "awk '{print $1" + '"-"' + "$2}' " + i + " > temporary1")

#                print("awk '{print $1" + '"-"' + "$4}' " + mapList + " | grep -n -x -w -F -f temporary1 | awk '{split($0, array, " + '":"' + "); print(array[1]) }' | awk 'FNR==NR{a[$1];next}FNR in a' - " + mapList + " > " + sample_id + ".map")
#                os.system("awk '{print $1" + '"-"' + "$4}' " + mapList + " | grep -n -x -w -F -f temporary1 | awk '{split($0, array, " + '":"' + "); print(array[1]) }' | awk 'FNR==NR{a[$1];next}FNR in a' - " + mapList + " > " + sample_id + ".map")
#                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "awk '{print $1" + '"-"' + "$4}' " + mapList + " | grep -n -x -w -F -f temporary1 | awk '{split($0, array, " + '":"' + "); print(array[1]) }' | awk 'FNR==NR{a[$1];next}FNR in a' - " + mapList + " > " + sample_id + ".map")
#                os.remove("temporary1")

#                pedout = open(sample_id + '.ped', 'w')
                csvLine2 = [sample_id]
                pileupin = open(i, 'r')     
                snpcount, discar = re.split(' ', os.popen('wc -l ' + i).read())
                snpcount = int(snpcount)
#                pedout.write(sample_id + ' ' + sample_id + ' 0 0 0 1 ')
#                pileupin = open(i, 'r')
                counter = 1
                countbads = 0
                for line in pileupin:
                    if counter <= snpcount:
                        line = line.rstrip()
                        if line.split("\t")[3] == "0":
                                csvLine2.append('NA')
                                countbads += 1
                        else:
                            chro, coord, refb, cov, rbase, qual = line.split('\t')
                            counter += 1
                            rbase = rbase.replace(".",refb.capitalize())
                            rbase = rbase.replace(",",refb.capitalize())
                            if int(cov) == 1:
                                if r"*" in rbase or r"+" in rbase or r"-" in rbase:
                                    csvLine2.append('NA')
                                elif r"$" in rbase:
                                    rbase = rbase[0]
                                    if exclTerminals == True: 
                                        csvLine2.append('NA')
                                    else:
                                        csvLine2.append(rbase.capitalize() + rbase.capitalize())

                                elif r"^" in rbase:
                                    rbase = rbase[2]
                                    if exclTerminals == True:
                                        csvLine2.append('NA')
                                    else:
                                        csvLine2.append(rbase.capitalize() +rbase.capitalize())
                                else:
                                    csvLine2.append(rbase.capitalize() + rbase.capitalize())
                            
                            ## If there is more than one call for this position
                            elif int(cov) > 1:
                                listBases = list(rbase)
                                ## Remove indels
                                while r"-" in listBases:
                                    for i in listBases:
                                        if i == r"-":
                                            indexi = listBases.index(r"-")
                                            delL = listBases[indexi+1]
                                            del listBases[indexi:(indexi+int(delL)+2)] # delete indel annotation
                                            del listBases[indexi-1] # delete indel
                                while r"+" in listBases:
                                    for i in listBases:
                                        if i == r"+":
                                            indexi = listBases.index(r"+")
                                            delL = listBases[indexi+1]
                                            del listBases[indexi:(indexi+int(delL)+2)] # delete indel annotation
                                            del listBases[indexi-1] # delete indel

                                ## Remove start and end of read calls is requested (^ and $)
                                while r"^" in listBases:
                                    for i in listBases:
                                        if i == r"^":
                                            indexi = listBases.index(r"^")
                                            del listBases[indexi] # delete ^
                                            del listBases[indexi] # delete the base quality after the ^
                                            if exclTerminals == True:
                                                del listBases[indexi] # delete the actual variant 
                                while r"$" in listBases:
                                    for i in listBases:
                                        if i == r"$":
                                            indexi = listBases.index(r"$")
                                            del listBases[indexi] # delete $
                                            if exclTerminals == True:
                                                del listBases[indexi-1] # delete the actual variant 

                                if len(listBases) != 0:
                                    chosenbase = random.choice(listBases)
                                    csvLine2.append(chosenbase.capitalize() + chosenbase.capitalize())
                                else:
                                    countbads += 1
                                    csvLine2.append('NA')
                    else:
                        pass

                writer.writerow(csvLine2)

    csvOut.close()

    print("Rscript " + hpsconvertonline_location + " " + cwd + "/" + 'HIrisPlex-S_toConvert.csv')
    os.system("Rscript " + hpsconvertonline_location + " " + cwd + "/" + 'HIrisPlex-S_toConvert.csv')

    #print("Rscript Testing.R " + cwd + "/" + 'HIrisPlex-S_toConvert.csv')
    #os.system("Rscript Testing.R " + cwd + "/" + 'HIrisPlex-S_toConvert.csv')












def CpGdamage(samIn="N18_merged.q30.rmdup.bam", cpgDatabase="GPL13534-11288_reduced_ultra_37.txt", CGasMissing=True,asFraction=True):
    """
    adasdasdasd
    """

    import pysam

    samfile = pysam.AlignmentFile(samIn, 'rb')
    cpgDB = open(cpgDatabase, 'r')
    cpgDict = {}

    cpgDB2 = cpgDB.readlines()[1:]
    print(len(cpgDB2))
    CpGmeth = {}
    hitCount = {}
    for ligna in cpgDB2:
        cgID, chro, pos, strand = ligna.split('\t')
        chrom = str("chr" + str(chro))
        strand = strand.rstrip('\n')
        cpgDict[chrom, pos] = cgID

        reads = samfile.fetch(chrom, int(pos) - 1, int(pos))  ##There seems to be an issue with the Illumina CpG database I downloaded, as coordinate of the C in a CpG seems to actually be -1 than shown
        for i in reads:
            # print(i)
            # print(int(pos)-1)    ##Position of desired CpG
            # print(str(i).split('\t')[3])    ##Starting position of read
            # print(str(i).split('\t')[9])    ##Read
            print("Methylated so far:", len(CpGmeth))
            print("Current chromosome: ", chrom, "\n")
            if int(pos) - 1 == int(str(i).split('\t')[3]):  ##Find if position of CpG is at the very start of the read
                if cgID in hitCount:  ##Counts how many 'CG's at start of reads exist for a specific cgID. Includes info for CpGs that either have damage pattern or nots. Needed this way for when CGasMissing == False
                    hitCount[cgID] += 1
                else:
                    hitCount[cgID] = 1
                ##Treats 'CG's at start of reads as missing information, because technically we can not be sure if they are methylated or not. They just are not damaged
                if CGasMissing == True:
                    if str(i).split('\t')[9][0] == "T":
                        # print("Methylated!")
                        if cgID in CpGmeth:
                            CpGmeth[cgID] += 1
                        else:
                            CpGmeth[cgID] = 1
                            # pass
                            # else:
                            #     print("Unmethylated!")


                            # print(CpGmeth)
                            # print(hitCount)
                            # print(hitCount[cgID])
                            # input("Got it")

                ##Assuming a 'CG' as unmethylated, which is technically not true since it simply does not have damage
                if CGasMissing == False:
                    CpGmeth[cgID] = 0  ##Statement defining line
                    if str(i).split('\t')[9][0] == "T":
                        # print("Methylated!")
                        if cgID in CpGmeth:
                            CpGmeth[cgID] += 1
                        else:
                            CpGmeth[cgID] = 1
                            # else:
                            #     print("Unmethylated!")

                            # print(CpGmeth)
                            # input("Got it")

    # print(CpGmeth)
    # print(hitCount)
    print("Len of hitcount: ", len(hitCount))
    print("Len of CpGmeth: ", len(CpGmeth))

    fileOut = open(str(samIn + ".betabinary"), 'w')
    if asFraction == True:
        fileOut.write("ID\tMethCount\tTotalHits\tBetaAsFrac\n")
        for item in CpGmeth:
            # print(item)
            # print(CpGmeth[item])
            beta = CpGmeth[item] / hitCount[item]
            fileOut.write(item + "\t" + str(CpGmeth[item]) + "\t" + str(hitCount[item]) + "\t" + str(beta) + "\n")
    elif asFraction == False:
        fileOut.write("ID\tMethCount\tTotalHits\tBeta\n")
        for item in CpGmeth:
            # print(item)
            # print(CpGmeth[item])
            if CpGmeth[item] > 1:
                beta = 1
            else:
                beta = CpGmeth[item]
            fileOut.write(item + "\t" + str(CpGmeth[item]) + "\t" + str(hitCount[item]) + "\t" + str(beta) + "\n")

    cpgDB.close()


def transFst(bimFile):
    """
    Outputs a list of SNPs characterized by having two transversion alleles (i.e. exclude CT and AG).

    :param bimFile:
    :return:
    """

    bimFileOpen = open(bimFile, "r")
    allFile = bimFileOpen.readlines()
    transList = []
    counter = 0
    for i in allFile:
        ct = str(i.split("\t")[4]) + str(i.split("\t")[5])
        if ct != "CT\n" and ct != "TC\n" and ct != "AG\n" and ct != "GA\n":
            transList.append(i.split("\t")[1])
        else:
            counter += 1

    with open(bimFile.split(".bim")[0] + "_transvOnly", "w") as fout:
        for i in transList:
            fout.write(i + "\n")


# print(len(transList))
#    print(counter)

def reorderBAM(index, extensionBAM=".q30.rmdup.bam", rename=True):
    """
    Reorders contigs in BAM files that do not match the order on the desired index.

    :param index:
    (str)   Path to index fasta file.
    :param extensionBAM:
    (str)   Extension of BAM files to find and work on. All files must have this extension. Only one string accepted.
    :param rename:
    (boolean)   Rename original BAM files by adding ".originalOrderBam", and then rename the newly ordered files to match the original file name. True or False.
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< REORDERING BAM CONTIGS - Reordering BAM file chromosomes according to given index genome >")
    fout.write("\n< REORDERING BAM CONTIGS - Reordering BAM file chromosomes according to given index genome >")

    ## Check if index genome already has a sequence dictionary
    counter = 0
    pathToIndex, indexFa = index.rsplit("/", 1)
    if re.search(".fasta$", indexFa) != None:
        indexNoFa = indexFa.rsplit(".fasta", 1)[0]
    elif re.search(".fa$", indexFa) != None:
        indexNoFa = indexFa.rsplit(".fa", 1)[0]
    for file in os.listdir(pathToIndex):
        if re.search(indexNoFa + ".dict", file) != None:
            counter += 1
    if counter == 0:
        print("\t> WARNING: Index dictionary file not found. Trying to generate it with PicardTools...")
        fout.write("\n\t> WARNING: Index dictionary file not found. Trying to generate it with PicardTools...")
        print("\t" + 'picard.jar CreateSequenceDictionary R=' + index + ' O=' + indexNoFa + '.dict')
        os.system('picard.jar CreateSequenceDictionary R=' + index + ' O=' + indexNoFa + '.dict')
        fout.write(time.strftime(
            "\n\t%H:%M:%S %Z ") + "\t" + 'picard.jar CreateSequenceDictionary R=' + index + ' O=' + indexNoFa + '.dict' + '\n')

    if counter == 0 and indexNoFa + ".dict" not in os.listdir(pathToIndex):
        print("\t> ERROR: Could not create or find index dictionary file. Aborting...")
        fout.write("\n\t> ERROR: Could not create or find index dictionary file. Aborting...")
        quit()

    ## Continue if all ok
    for file in os.listdir(cwd):
        if re.search(extensionBAM, file) != None and re.search(".bam$", file) != None:
            print("\t" + 'picard.jar ReorderSam R=' + index + ' O=' + file + '.reorder I=' + file)
            os.system('picard.jar ReorderSam R=' + index + ' O=' + file + '.reorder I=' + file)
            fout.write(time.strftime(
                "\n\t%H:%M:%S %Z ") + "\t" + 'picard.jar ReorderSam R=' + index + ' O=' + file + '.reorder I=' + file + '\n')
            ## Rename if True
            if rename == True:
                os.rename(file, file + ".originalOrderBam")
                os.rename(file + ".reorder", file)


def editHeaderBAMbatch(extensionBAM=".q30.rmdup.bam", RGid="1", RGlb="lib1", RGpl="ILLUMINA", RGpu="unit1",rename=True):
    """
    Edits/adds headers for/to all BAM files in folder with given extension.

    :param extensionBAM:
    (str)   Extension of BAM files to find and work on. All files must have this extension. Only one string accepted.
    :param RGid:
    (str)   Group tag. Default is "1"
    :param RGlb:
    (str)   Library tag. Default is "lib1"
    :param RGpl:
    (str)   Sequencing platform tag. Default is "ILLUMINA"
    :param RGpu:
    (str)   Platform unit tag. Default is "unit1"
    :param rename:
    (boolean)   Rename original BAM files by adding ".originalOrderBam", and then rename the newly ordered files to match the original file name. True or False.
    :return:
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< EDITING BAM HEADERS - Edits BAM file headers >")
    fout.write("\n< EDITING BAM HEADERS - Edits BAM file headers >")

    for file in os.listdir(cwd):
        if re.search(extensionBAM, file) != None and re.search(".bam$", file) != None:
            sample_id, sep = file.split(extensionBAM, 1)
            print(sample_id)
            print(
                "\t" + 'picard.jar AddOrReplaceReadGroups I=' + file + ' O=' + file + '.newheader RGID=' + RGid + ' RGLB=' + RGlb + ' RGPL=' + RGpl + ' RGPU=' + RGpu + ' RGSM=' + sample_id)
            os.system(
                'picard.jar AddOrReplaceReadGroups I=' + file + ' O=' + file + '.newheader RGID=' + RGid + ' RGLB=' + RGlb + ' RGPL=' + RGpl + ' RGPU=' + RGpu + ' RGSM=' + sample_id)
            fout.write(time.strftime(
                "\n\t%H:%M:%S %Z ") + "\t" + 'picard.jar AddOrReplaceReadGroups I=' + file + ' O=' + file + '.newheader RGID=' + RGid + ' RGLB=' + RGlb + ' RGPL=' + RGpl + ' RGPU=' + RGpu + ' RGSM=' + sample_id + '\n')
            ## Rename if True
            if rename == True:
                os.rename(file, file + ".originalHeaderBam")
                os.rename(file + ".newheader", file)


def test_su():
    import sys
    prac3functions = sys.argv[0]
    a, b, c, d = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    mylist = [a, b, c, d]
    print(mylist)

    def calculate_sum(a):
        summ = sum(a)
        print("The sum of the four integers {} is: {}".format(a, summ))

    def calculate_product(b):
        product = ((b[0]) * (b[1]) * (b[2]) * (b[3]))
        print("The product of the four integers {} is: {}".format(b, product))

    calculate_product(mylist)
    calculate_sum(mylist)


def angsd_contamination(extensionBAM,mode="C"):
    """

    :param extensionBAM:

    :param mode:
    Use R or C mode of contamiantion estimation. Default C.
    """

    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< ANGSD - Estimating contamination on male individuals with ANGSD >")
    fout.write("\n< ANGSD - Estimating contamination on male individuals with ANGSD >")

    for file in os.listdir(cwd):
        if re.search(extensionBAM, file) != None and re.search(".bam$", file) != None:
            sample_id, sep = file.split(extensionBAM, 1)
            print("\t/home/dfernandes/Software/angsd/angsd -i "+file+" -r chrX:5000000-154900000 -doCounts 1  -iCounts 1 -minMapQ 30 -minQ 20")
            os.system("/home/dfernandes/Software/angsd/angsd -i "+file+" -r chrX:5000000-154900000 -doCounts 1  -iCounts 1 -minMapQ 30 -minQ 20")
            os.rename("angsdput.icnts.gz",sample_id+".icnts.gz")
            if mode == "C":
                print("\tnohup /home/dfernandes/Software/angsd/misc/contamination -a "+sample_id+".icnts.gz -h /home/dfernandes/Software/angsd/RES/HapMapChrX.gz > "+sample_id+"angsd")
                os.system("nohup /home/dfernandes/Software/angsd/misc/contamination -a "+sample_id+".icnts.gz -h /home/dfernandes/Software/angsd/RES/HapMapChrX.gz > "+sample_id+"angsd")
            elif mode == "R":
                print("nohup Rscript /home/dfernandes/Software/angsd/R/contamination.R mapFile=/home/dfernandes/Software/angsd/RES/chrX.unique.gz hapFile=/home/dfernandes/Software/angsd/RES/HapMapChrX.gz countFile=" + sample_id + ".icnts.gz > " + sample_id + "angsd")
                os.system("nohup Rscript /home/dfernandes/Software/angsd/R/contamination.R mapFile=/home/dfernandes/Software/angsd/RES/chrX.unique.gz hapFile=/home/dfernandes/Software/angsd/RES/HapMapChrX.gz countFile=" + sample_id + ".icnts.gz > " + sample_id + "angsd")

def damageRestrict(samIn, refGenome):
    """
    adasdasdasd
    """

    import subprocess
    import pysam

    # import time
    # start_time = time.time()

    samfile = pysam.AlignmentFile(samIn, 'rb')
    samOut = pysam.AlignmentFile(samIn+".damageRestrict", "wb", template=samfile)

    ## Create chr order dictionary
    os.system("samtools idxstats "+samIn+" > temp1")
    temp1open = open("temp1", "r")
    chrOrder = {}
    counter = 0
    for line in temp1open:
        chrOrder[counter] = line.split()[0]
        counter += 1
    os.remove("temp1")

    lista = []
    counter = 1
    for i in samfile:
        if str(i).split()[9][0] == "T":
            lista.append((counter))
        counter += 1
    samfile.close()
    samfile = pysam.AlignmentFile(samIn, 'rb')

    #countCs = 0 ## If want to get a percentage, like mapDamage
    counter = 1
    for i in samfile:
        if counter in lista:
            chro = int(str(i).split()[2])
            pos = int(str(i).split()[3])+1  ## Pysam reports 1bp less on the coordinates, so adding +1
            readBp1 = str(i).split()[9][0]

            proc = subprocess.Popen(["samtools faidx "+refGenome+" "+chrOrder[chro]+":"+str(pos)+"-"+str(pos)],stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
            refBase = str(out.split()[1])[2].upper()

            # if readBp1 == "C":  ## If want to get a percentage, like mapDamage
            #     countCs += 1

            if readBp1 == "T" and refBase == "C":
                samOut.write(i)
            counter +=1
        else:
            counter +=1


    samfile.close()
    samOut.close()

    #print("--- %s seconds ---" % (time.time() - start_time))
    #print(countCs)

def mergeBamsBatch(qualThreshold = 30, remSourceFiles = False):
    """
    Merges each sample's 'rmdup.bam' files originated by multiple lanes.

    Input:
    Foo1_L001.q30.rmdup.bam
    Foo1_L002.q30.rmdup.bam
    Foo1_L003.q30.rmdup.bam
    Foo1_L004.q30.rmdup.bam

    Output:
    Foo1_merged.q30.rmdup.bam
    Foo1_merged.q30.rmdup.bai

    :param remSourceFiles:
    (boolean)   Remove source files. Default = False.
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< MERGEBAMSBATCH - Batch merging each sample's bam files >")
    fout.write("\n< MERGEBAMSBATCH - Batch merging each sample's bam files >")

    qualThreshold = str(qualThreshold)
    bamDict = {}
    for i in os.listdir(cwd):
        #print(i)
        ja = re.search(r'_S\d\d_L\d\d\d\.q..\.rmdup\.bam$', i)
        ba = re.search(r'\.bai$', i)
        if ja is not None and ba is None:
            sampleName = i[:ja.start()]
            #print(sampleName)

            if sampleName in bamDict.keys():
                if i in bamDict[sampleName]:
                    pass

                else:
                    bamDict[sampleName].append(i)
                    #print("b")
            else:
                bamDict[sampleName] = [i]

    for i in os.listdir(cwd):
        # print(i)
        ja = re.search(r'_S\d_L\d\d\d\.q..\.rmdup\.bam$', i)
        ba = re.search(r'\.bai$', i)
        if ja is not None and ba is None:
            sampleName = i[:ja.start()]
            # print(sampleName)

            if sampleName in bamDict.keys():
                if i in bamDict[sampleName]:
                    pass

                else:
                    bamDict[sampleName].append(i)
                    # print("b")
            else:
                bamDict[sampleName] = [i]

            #print(ja.group())
            #print(ja.start())
            #print("\n")
    print(bamDict)
    print(bamDict.keys())

    for i in bamDict.keys():
        #print(i)
        command = "picard.jar MergeSamFiles "
        print("\t> Merging all bam files for sample with id " + i + " and sorting")
        fout.write("\n\n\t> Merging all bam files for sample with id " + i + " and sorting")
        for a in bamDict[i]:
            #print(a)
            command = command + "I=" + a + " "
        command = command + "O=" + i + "_merged.q" + qualThreshold + ".bam USE_THREADING=True"
        print(command)
        os.system(command)
        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + command)

        ## Move merged files into new folder
        if "mergedlanes" not in os.listdir(cwd):
            os.system('mkdir -p '+cwd+'/mergedlanes')
        os.system('mv ' + i + '_merged.q' + qualThreshold + '.bam ' + cwd + '/' + "mergedlanes/")

        os.chdir(cwd + '/mergedlanes')

        print("\tsamtools rmdup -s " + i + "_merged.q" + qualThreshold + ".bam" + " " + i + "_merged.q" + qualThreshold + ".rmdup.bam")
        os.system("samtools rmdup -s " + i + "_merged.q" + qualThreshold + ".bam" + " " + i + "_merged.q" + qualThreshold + ".rmdup.bam")
        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools rmdup -s " + i + "_merged.q" + qualThreshold + ".bam" + " " + i + "_merged.q" + qualThreshold + ".rmdup.bam")

        print("\tsamtools flagstat " + i + "_merged.q" + qualThreshold + ".rmdup.bam" + " > " + i + "_merged.q" + qualThreshold + ".summary.txt")
        os.system("samtools flagstat " + i + "_merged.q" + qualThreshold + ".rmdup.bam" + " > " + i + "_merged.q" + qualThreshold + ".summary.txt")
        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools flagstat " + i + "_merged.q" + qualThreshold + ".rmdup.bam" + " > " + i + "_merged.q" + qualThreshold + ".summary.txt")
        print("\tsamtools index " + i + "_merged.q" + qualThreshold + ".rmdup.bam" + " " + i + "_merged.q" + qualThreshold + ".rmdup.bam.bai")
        os.system("samtools index " + i + "_merged.q" + qualThreshold + ".rmdup.bam" + " > " + i + "_merged.q" + qualThreshold + ".rmdup.bam.bai")
        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools index " + i + "_merged.q" + qualThreshold + ".rmdup.bam" + " > " + i + "_merged.q" + qualThreshold + ".rmdup.bam.bai")
        print("\tsamtools idxstats " + i + "_merged.q" + qualThreshold + ".rmdup.bam" + " > " + i + "_merged.chr.txt")
        os.system("samtools idxstats " + i + "_merged.q" + qualThreshold + ".rmdup.bam" + " > " + i + "_merged.chr.txt")
        fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools idxstats " + i + "_merged.q" + qualThreshold + ".rmdup.bam" + " > " + i + "_merged.chr.txt")

        os.chdir(cwd)

def calico(refLoc,calLoc,qualThreshold = 30,multipleLanes = False):
    """
    Estimates mitochondrial contamination in human samples from "trimmed.fastq" files, using the Calico script by Pontus Skoglund.

    Input:
    Foo1_S01_L001.trimmed.fastq
    Foo1_S01_L002.trimmed.fastq
    Foo2_S02_L001.trimmed.fastq
    Foo2_S02_L001.trimmed.fastq

    Output:
    Foo1_S01.trimmed.fastq
    Foo2_S02.trimmed.fastq
    Foo1_S01_Calico/.q30.bam/.q30.rmdup.bam/.q30.rmdup.bam.bai/.q30.sort.bam/.sai/.sam
    Foo2_S02_Calico/.q30.bam/.q30.rmdup.bam/.q30.rmdup.bam.bai/.q30.sort.bam/.sai/.sam
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< CALICO - Estimating mitochondrial contamination >")
    fout.write("\n< CALICO - Estimating mitochondrial contamination >")

    qualThreshold = str(qualThreshold)

    if seqtkSim == 1:
        fileExp = r'seqtk\.fastq$'
    else:
        fileExp = r'trimmed\.fastq$'

    bamDict = {}
    for i in os.listdir(cwd):
        ja = re.search(fileExp, i)
        #ba = re.search(r'\.bai$', i)
        if ja is not None:
            sampleName = i[:ja.start()]
            sample_id = str.split(sampleName,"%")[0]
            #sample_id = sample_id[:-5] ## Possible fix, need to confirm

            if sample_id in bamDict.keys():
                if i in bamDict[sample_id]:
                    pass

                else:
                    bamDict[sample_id].append(i)
            else:
                bamDict[sample_id] = [i]

    print(bamDict)

    for it in bamDict:
        print(it)
        if multipleLanes == True:
            commandCat = "cat "
            for i in bamDict[it]:
                print(i)
                commandCat += str(i) + (" ")
            commandCat += " > " + it + ".trimmed.fastq"

            i = it + ".trimmed.fastq"
            #sample_id = it
            print(commandCat)
            os.system(commandCat)
        #print(i)
        if multipleLanes == True:
            i = it + ".trimmed.fastq"
        elif multipleLanes == False:
            i = bamDict[it][0]
        sample_id = it
        print("\t> Aligning trimmed sequences from " + "'" + i + "'" + " to reference genome...")
        print('\tbwa aln -l 1000 -t 2 ' + refLoc + ' ' + i + ' > ' + sample_id + '_Calico.sai' + ' 2> ' + sample_id + '_Calico.bwa_aln.log' + "\n")
        os.system('bwa aln -l 1000 -t 2 ' + refLoc + ' ' + i + ' > ' + sample_id + '_Calico.sai' + ' 2> ' + sample_id + '_Calico.bwa_aln.log')
        fout.write("\n\t> Aligning trimmed sequences from " + "'" + i + "'" + " to reference genome..." + '\n\tbwa aln -l 1000 -t 2 ' + refLoc + ' ' + i + ' > ' + sample_id + '_Calico.sai' + ' 2> ' + sample_id + '_Calico.bwa_aln.log' + "\n")

        it=i
        i= sample_id + "_Calico.sai"
        print("\t> Converting " + i + " to '.sam' and adding tags...")
        print("\t" + r"bwa samse -r '@RG\tID:" + "calico" + r"\tSM:" + sample_id + r"\tPL:ILLUMINA' " + refLoc + ' ' + i + ' ' + it + ' > ' + sample_id + '_Calico.sam')
        os.system(r"bwa samse -r '@RG\tID:" + "calico" + r"\tSM:" + sample_id + r"\tPL:ILLUMINA' " + refLoc + ' ' + i + ' ' + it + ' > ' + sample_id + '_Calico.sam')
        fout.write("\n\t> Converting " + i + " to '_Calico.sam' and adding tags..." + time.strftime("\n\t%H:%M:%S %Z ") + "\t" + r"bwa samse -r '@RG\tID:" + "calico" + r"\tSM:" + sample_id + r"\tPL:ILLUMINA' " + refLoc + ' ' + i + ' ' + it + ' > ' + sample_id + '_Calico.sam' + '\n')
        print("\t> Converting " + sample_id + '_Calico.sam' + " to '.q" + qualThreshold + ".bam'...")
        print("\t"r"samtools view -Sb -q " + qualThreshold + " -F 4 " + sample_id + "_Calico.sam" + " > " + sample_id + '_Calico.q' + qualThreshold + '.bam')
        os.system(r"samtools view -Sb -q " + qualThreshold + " -F 4 " + sample_id + "_Calico.sam" + " > " + sample_id + '_Calico.q' + qualThreshold + '.bam')
        fout.write("\t> Converting " + sample_id + '_Calico.sam' + " to '_Calico.q' + qualThreshold + '.bam'..." + time.strftime("\n\t%H:%M:%S %Z ") + "\t"r"samtools view -Sb -q " + qualThreshold + " -F 4 " + sample_id + "_Calico.sam" + " > " + sample_id + '_Calico.q' + qualThreshold + '.bam' + '\n')
        print("\t")



        print("\t" + 'samtools sort ' + sample_id + '_Calico.q' + qualThreshold + '.bam' + ' -o ' + sample_id + '_Calico.q' + qualThreshold + '.sort.bam' + "\n")
        os.system('samtools sort ' + sample_id + '_Calico.q' + qualThreshold + '.bam' + ' -o ' + sample_id + '_Calico.q' + qualThreshold + '.sort.bam')
        fout.write("\n\t> Sorting " + sample_id + '_Calico.q' + qualThreshold + '.bam' + "..." + time.strftime("\n\t%H:%M:%S %Z ") + '\tsamtools sort ' + i + ' -o ' + sample_id + '_Calico.q' + qualThreshold + '.sort.bam' + "\n")

        print("\t" + 'samtools rmdup -s ' + sample_id + '_Calico.q' + qualThreshold + '.sort.bam' + ' ' + sample_id + '_Calico.q' + qualThreshold + '.rmdup.bam')
        os.system('samtools rmdup -s ' + sample_id + '_Calico.q' + qualThreshold + '.sort.bam' + ' ' + sample_id + '_Calico.q' + qualThreshold + '.rmdup.bam')
        fout.write("\n\t> Removing duplicates from file " + i + "..." + time.strftime("\n\t%H:%M:%S %Z ") + "\t" + 'samtools rmdup -s ' + i + ' ' + sample_id + '.rmdup.bam' + '\n')

        print("\t> Indexing '*.bam' file {}...".format(sample_id))
        print('samtools index ' + sample_id + '_Calico.q' + qualThreshold + '.rmdup.bam' + ' ' + sample_id + '_Calico.q' + qualThreshold + '.rmdup.bam' + '.bai')
        os.system('samtools index ' + sample_id + '_Calico.q' + qualThreshold + '.rmdup.bam' + ' ' + sample_id + '_Calico.q' + qualThreshold + '.rmdup.bam' + '.bai')
        print("\t")

        print('samtools mpileup -q' + qualThreshold +' -Q' + qualThreshold + ' -f ' + refLoc + ' ' + sample_id + '_Calico.q' + qualThreshold + '.rmdup.bam | python ' + calLoc + ' --maxdepth 100000 --indels > ' + sample_id + '_Calico')
        os.system('samtools mpileup -q' + qualThreshold +' -Q' + qualThreshold + ' -f ' + refLoc + ' ' + sample_id + '_Calico.q' + qualThreshold + '.rmdup.bam | python ' + calLoc + ' --maxdepth 100000 --indels > ' + sample_id + '_Calico')

    ## Collect and compile data
    calOut = open("CalicoSummary.txt","w")
    calOut.write("SAMPLE\tSITES\tMAJ\tMIN\tCONT\tCI_low\tCI_up\n")
    for i in os.listdir(cwd):
        if re.search(r'_Calico$', i) is not None:
            print(i)

            calIn = open(i,'r')
            calIn.readline()
            calOut.write(i.strip(r'_Calico$'))
            calOut.write("\t")
            calOut.write(calIn.readline())
            #print(calIn.readline())
    calOut.close()

def fastqdemux(inFastq,sampleList,mismatch = 0):
    """
    Demultiplexes unindexed Illumina reads, such as the ones stored in /Projects/Unindexed Reads/Samples/ in the user's Basespace account.
    Always provide the reverse complement of the forward adapter sequences for the I7 and I5.
    Paired-end sequencing not supported yet.

    Current version average writing speed:
    130MB/minute - 2.22MB/second

    :param inFastq:
    FASTQ(.GZ) file to parse. This script was prepared to work with "Undetermined" FASTQ files from Basespace but can accept any FASTQ files with the index sequences at the end of the first line. Example:
    @NB501111:90:HYN5FBGX5:1:11101:18918:1050 1:N:0:CTCGATG+CTCGATG (for dual-indexed libraries)
    @NB501111:90:HYN5FBGX5:1:11101:18918:1050 1:N:0:CTCGATG (for single-indexed libraries)
    The filename must contain lane number information ("L00x"). Example:
    Undetermined_S0_L004_R1_001.fastq(.gz)

    :param sampleList:
    A tab-spaced text file with sample name(s), reverse complement of index1 (I7) sequence, and reverse complement of index2 (I5) sequence (when dual-indexed), respectively in columns one, two, and three. Example for dual-indexed libraries:
    Sample1\tAGTCGTT\tTGCCTAT
    Sample2\tGCTTCTA\tCGGGTAC

    :param mismatch:
    Number of allowed mismatches (default = 0)
    """
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< DEMUXQ - Demultiplexing FASTQ into individual sample files >")
    fout.write("\n< DEMUXQ - Demultiplexing FASTQ into individual sample files >")

    sampIn = open(sampleList, 'r')
    sampDict = {}


    if(len(sampIn.readline().split("\t"))) == 2:
        sampIn.seek(0)

        for line in sampIn:
            if str.split(line)[0] in sampDict.values():
                print("\t> ERROR: Found duplicated index ",str.split(line)[0])
                print("\t> Please ensure all index sequences are unique. Quitting now...")
                quit()

            if str.split(line)[1] in sampDict.keys():
                print("\t> ERROR: Found duplicated sample ",str.split(line)[1])
                print("\t> Please ensure all sample names are unique. Quitting now...")
                quit()
            else:
                sampDict[str.split(line)[1]] = str.split(line)[0]

    elif(len(sampIn.readline().split("\t"))) == 3:
        sampIn.seek(0)

        for line in sampIn:
            if str.split(line)[0] in sampDict.values():
                print("\t> ERROR: Found duplicated index ",str.split(line)[0])
                print("\t> Please ensure all index sequences are unique. Quitting now...")
                quit()

            if str.split(line)[1] in sampDict.keys():
                print("\t> ERROR: Found duplicated sample ",str.split(line)[1])
                print("\t> Please ensure all sample names are unique. Quitting now...")
                quit()
            else:
                sampDict[str.split(line)[1]+"\\\+"+str.split(line)[2]] = str.split(line)[0]

    sampIn.close()
    sampList=list(sampDict.keys())

    ## Check for previous files and delete them
    lnum = inFastq.split("L00")[1][0]
    for i in sampList:
        for items in os.listdir(cwd):
            if sampDict[i] in items:
                try:
                    os.remove(items)
                except:
                    pass
    try:
        os.remove("Unindexed_S0_" + "L00" + lnum + "_R1_001.fastq")
    except:
        pass

    ## Create results-related dicts
    sampNum = {}
    seqa = 1
    for i in sampDict.keys():
        sampNum[i] = seqa
        seqa += 1

    sampFinds = {}
    for i in sampDict.keys():
        sampFinds[i] = 0
    countUndet = 0

    ## Iterate over the FASTQ(.GZ) file reads
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    if inFastq[-3:] == "stq":
        with open(inFastq) as in_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                ite = 0
                found = title.split(':')[-1].replace("+", "\\\+")

                ## The next loop is the one that slows down the whole process

                for onekey in sampDict.keys():
                    if len([b for a, b in zip(onekey, found) if b != a]) <= mismatch:  # shows mismatched bases from Found
                    #if len(regex.findall("(" + found + "){e<=" + str(mismatchCompensate) + "}",onekey)) > 0:
                        output_file = (sampDict[onekey] + "_S" + str(sampNum[onekey]) + "_L00" + lnum + "_R1_001.fastq")
                        with open(output_file, "a") as out_handle:
                            out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                            sampFinds[onekey] += 1
                            ite = 1
                if ite == 0:
                    output_file = ("Unindexed_S0_" + "L00" + lnum + "_R1_001.fastq")
                    with open(output_file, "a") as out_handle:
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        countUndet += 1

    if inFastq[-3:] == ".gz":
        with gzip.open(inFastq, "rt") as in_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                ite = 0
                found = title.split(':')[-1].replace("+", "\\\+")
                for onekey in sampDict.keys():
                    if len([b for a, b in zip(onekey, found) if b != a]) <= mismatch:  # shows mismatched bases from Found
                    #if len(regex.findall("(" + found + "){e<=" + str(mismatchCompensate) + "}",onekey)) > 0:
                        output_file = (sampDict[onekey] + "_S" + str(sampNum[onekey]) + "_L00" + lnum + "_R1_001.fastq")
                        with open(output_file, "a") as out_handle:
                            out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                            sampFinds[onekey] += 1
                            ite = 1
                if ite == 0:
                    output_file = ("Unindexed_S0_" + "L00" + lnum + "_R1_001.fastq")
                    with open(output_file, "a") as out_handle:
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                        countUndet += 1

    ## Create empty FASTQ files for unlocated samples
    import glob
    for i in list(sampDict.values()):
        if len(glob.glob(i+"*")) == 0:
            with open((i + "_S" + str(list(sampDict.values()).index(i) + 1) + "_L00" + lnum + "_R1_001.fastq"), "w"):
                pass

    ## Print short report
    print("\tHits found with max "+str(mismatch)+" mismatch(es):")
    for i in sampFinds:
        print("\t\t"+sampDict[i]+" - "+str(sampFinds[i]))
    print("\n\t\tUndetermined/no matches - "+str(countUndet))


def mergedEnd_Analysis():
    """
    Retrieves data from NextSeq's merged lanes, per sample. Requires "Reads_Summary.txt" file produced by function "end_analysis()".

    Output:
    ReadsMergeSummary.txt
    """

    os.chdir(originalWD)
    cwd = os.getcwd()
    script = open("RmergeLaneData.R", 'w')
    script.write('setwd("')
    script.write(cwd)
    script.write('")')

    script.write('\nwdfiles = list.files(getwd())')
    script.write('\ntbl = read.table("Reads_Summary.txt", sep="\\t", fill=T, header=T, stringsAsFactors=F)')
    script.write('\nsample_list = c()\nfor(i in list.files(pattern="_S[:0-9:]{1,2}_L[:0-9:]{3}.q[0-9][0-9].summary.txt$")) {\n\tsample_list = c(sample_list, strsplit(i,"_S[:0-9:]{1,2}_L[:0-9:]{3}.q[0-9][0-9].summary.txt$")[[1]][1])\n}\nsample_list = unique(sample_list)')
    script.write('\ntblOut = data.frame(matrix(nrow = length(sample_list),ncol = length(colnames(tbl))))\ncolnames(tblOut) = colnames(tbl)\ntmpTbl = data.frame(matrix(nrow = 4,ncol = length(colnames(tbl))))')
    script.write('\nitSamp = 1\n\nfor(i in sample_list){\n\titLane = 1\n\tfor(lin in grep(paste0("^",i,"_S[:0-9:]{1,2}_L[:0-9:]{3}$"), tbl$Sample)) {\n\t\ttmpTbl[itLane,] = tbl[lin,]\n\t\titLane = itLane +1\n\t}\n\ttblOut[itSamp,1] = i\n\ttblOut[itSamp,2] = sum(tmpTbl$X2)\n\ttblOut[itSamp,3] = median(tmpTbl$X3)\n\ttblOut[itSamp,4] = sum(tmpTbl$X4)\n\ttblOut[itSamp,5] = sum(tmpTbl$X5)\n\titSamp = itSamp +1\n}')
    script.write('\ntblOut[,6] = format(round((tblOut$Aligned_reads/tblOut$Total_reads)*100,2))\nmergedFiles = list.files(path = paste0(getwd(),"/mergedlanes"),pattern = "\\\merged.q[0-9][0-9].summary.txt$")')
    script.write('\nfor(i in mergedFiles) {\n\tsmpName = as.character(strsplit(i,"_merged.q[0-9][0-9].summary.txt"))\n\tdamageP = NA\n\ttblOut[grep(paste0("^",smpName,"$"),tblOut$Sample),7] = strsplit(read.table(paste0(getwd(),"/mergedlanes/",i),sep = "\\t",fill = T,header = F,stringsAsFactors = F)[1,]," ")[[1]][1]\n\tCTvalue = NA\n\tfiveCT = read.table(paste0(getwd(),"/mergedlanes/mapdamage/mapdamage_",smpName,"_merged/5pCtoT_freq.txt"),header = T,stringsAsFactors = F)[[2]]\n\tif(fiveCT[1] > fiveCT[2] && fiveCT[2] > fiveCT[3] && fiveCT[1] > fiveCT[length(fiveCT)]){\n\t\tCTvalue = fiveCT[1]\n\t}')
    script.write('\n\tGAvalue = NA\n\tthreeGA = read.table(paste0(getwd(),"/mergedlanes/mapdamage/mapdamage_",smpName,"_merged/3pGtoA_freq.txt"),header = T,stringsAsFactors = F)[[2]]\n\tif(threeGA[1] > threeGA[2] && threeGA[2] > threeGA[3] && threeGA[1] > threeGA[length(threeGA)]){\n\t\tGAvalue = threeGA[1]\n\t}')
    script.write('\n\tif(is.na(CTvalue) == F && is.na(GAvalue) == F){\n\t\tdamageP = paste0(format(round(CTvalue,3),nsmall=3)," | ",format(round(GAvalue,3),nsmall=3))\n\t}\n\ttblOut[grep(paste0("^",smpName,"$"),tblOut$Sample),10] = damageP')
    script.write('\n\tsex = try(read.table(paste0(getwd(),"/mergedlanes/sexing_",smpName,"_merged.txt"),header = T,sep="\\t",stringsAsFactors = F),silent=T)\n\tif(is.data.frame(sex)){\n\t\ttblOut[grep(paste0("^",smpName,"$"),tblOut$Sample),11] = sex$Assignment\n\t}\n}\n\ntblOut$X.Endogenous_Unique = format(round((as.numeric(tblOut$Unique_Reads)/tblOut$Total_reads)*100,2))\n\nwrite.table(x=tblOut,file="ReadsMergeSummary.txt",quote = F,sep = "\t",col.names = T,row.names = F)')
    script.close()

    #os.system("R CMD BATCH RmergeLaneData.R &")


def concatenateFQGZ():
    """
    Concatenates lane data from a Next/HiSeq run. Sample id is defined as the file name before "_L00X_R1.fastq.gz".

    Output:
    ReadsMergeSummary.txt
    """

    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< CONCATENATEFQGZ - Concatenating *fastq.gz* lane data from a Next/HiSeq run by sample >")
    fout.write("\n< CONCATENATEFQGZ - Concatenating *fastq.gz* lane data from a Next/HiSeq run by sample >")

    sampleDict = {}
    for i in os.listdir(cwd):
        ja = re.search(r'_S\d\d_L\d\d\d_R\d_00\d\.fastq\.gz$', i)
        if ja is not None:
            sample_id = i[:ja.start()]
            if sample_id in sampleDict.keys():
                if i in sampleDict[sample_id]:
                    pass
                else:
                    sampleDict[sample_id].append(i)
            else:
                sampleDict[sample_id] = [i]

    for it in sampleDict:
        commandCat = "cat "
        for i in sampleDict[it]:
            commandCat += str(i) + (" ")
        commandCat += " > " + it + "conc_S01_L001_R1_001.fastq.gz"

        print(commandCat)
        os.system(commandCat)
        fout.write("\t"+commandCat)



def blastnQuery(queryLink):
    """



    :param queryLink:
    :return:
    """

    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< BLASTNQUERY - Batch alignment with BLAST nucleotide query functionality >")
    fout.write("\n< BLASTNQUERY - Batch alignment with BLAST nucleotide query functionality >")





def imputeWithBeagle(curateRefPanel = False, bcftoolsNthreads = 1):



    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< IMPUTAION W BEAGLE - Impute genotypes using Beagle >")
    fout.write("\n< IMPUTAION W BEAGLE - Impute genotypes using Beagle >\n")



    if curateRefPanel == True:
        for i in os.listdir(cwd):
            if re.search('\.vcf\.gz$', i) is not None:
                #print(i)
                sample_id, sep = i.split(".vcf.gz")
                #print(sample_id)

                ## Keeping only biallelic SNPs
                print('\t> Keeping only bi-allelic SNPs for file: ', i)
                fout.write('\t> Keeping only bi-allelic SNPs for file: ' + i + "\n")
                print('bcftools view -m2 -M2 -v snps '+ i + ' --output-file ' +sample_id + '_bialSNPs.vcf.gz --output-type z --threads ' + str(bcftoolsNthreads))
                #os.system('bcftools view -m2 -M2 -v snps '+ i + ' --output-file ' +sample_id + '_bialSNPs.vcf.gz --output-type z --threads ' + str(bcftoolsNthreads))
                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + 'bcftools view -m2 -M2 -v snps '+ i + ' --output-file ' +sample_id + '_bialSNPs.vcf.gz --output-type z --threads ' + str(bcftoolsNthreads))

        ## Move old files into separate folder
            ## Checking if folder base_files exists, otherwise, create it
        newdir = cwd + '/' + "base_files"
        if 'base_files' in os.listdir(cwd):
            pass
        else:
            os.mkdir(newdir)

        for i in os.listdir(cwd):
            if re.search('.vcf.gz$', i) != None and '_bialSNPs' not in i:
                os.system('mv ' + i + ' ' + newdir)
                os.system('mv ' + i + '.tbi ' + newdir)

        ##### STATS
        indSNPs = []
        summaryOut = open('Summary_ALL', 'w')
        summaryOut.write("Sample\tNumber of SNPs\n")
        ## Run and compile stats
        for i in os.listdir(cwd):
            if re.search('_bialSNPs.vcf.gz$', i) != None:
                print('\t> Running bcftools stats for file: ', i)
                fout.write('\t> Running bcftools stats for file: ' + i + "\n")
                sample_id_stats, sep = i.split(".vcf.gz")
                print('bcftools stats ' + i + ' > Summary_stats_' + sample_id_stats)
                #os.system('bcftools stats ' + i + ' > Summary_stats_' + sample_id_stats)
                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + 'bcftools stats ' + i + ' > Summary_stats_' + sample_id_stats)
        ## Compile stats from all files
        print('\t> Compiling numbers of SNPs from summary data')
        fout.write('\t> Compiling numbers of SNPs from summary data')
        for i in os.listdir(cwd):
            if re.search('^Summary_stats_', i) != None:
                sep, sample_id_summary = i.split("Summary_stats_")
                readSummary = open(i,'r').readlines()
                numSnps = readSummary[[readSummary.index(s) for s in readSummary if "\tnumber of SNPs:\t" in s][0]].split("\t")[-1]
                summaryOut.write(sample_id_summary + "\t" + numSnps)
                indSNPs.append(int(numSnps))

        ## Create SNP position list
        print('\t> Creating SNP interval list file')
        fout.write('\t> Creating SNP interval list file')
        intervalFile = open('ALL_intervals.list','w')
        for i in os.listdir(cwd):
            if re.search('_bialSNPs.vcf.gz$', i) != None:
                sample_id_stats, sep = i.split("_bialSNPs.vcf.gz")
                print('zcat ' + i + r" | grep -v \# | awk '" + '{print $1":"$2"-"$2}' + "' >> ALL_intervals.list")
                os.system('zcat ' + i + r" | grep -v \# | awk '" + '{print $1":"$2"-"$2}' + "' >> ALL_intervals.list")
                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + 'zcat ' + i + r" | grep -v \# | awk '" + '{print $1":"$2"-"$2}' + "' >> ALL_intervals.list")
        intervalFile.close()

    totalSNPs = sum(indSNPs)
    summaryOut.write("TOTAL_SNPs\t" + str(totalSNPs))
    summaryOut.close()


def batchBAMtoFASTQ(picard = False):
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< BATCH BAM TO FASTQ - Convert files >")
    fout.write("\n< BATCH BAM TO FASTQ - Convert files >\n")

    for i in os.listdir(cwd):
        if re.search('.bam$', i) != None:
            if picard == False:
                print('\t> Running bedtools converter for file: ', i)
                fout.write('\t> Running bedtools converter for file: ' + i + "\n")
                sample_id, sep = i.split(".bam")
                #print(sample_id)
                print("bedtools bamtofastq -i " + i + " -fq " + sample_id + ".fastq")
                os.system("bedtools bamtofastq -i " + i + " -fq " + sample_id + ".fastq")
                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "bedtools bamtofastq -i " + i + " -fq " + sample_id + ".fastq")
            if picard == True:
                print('\t> Running PicardTools converter for file: ', i)
                fout.write('\t> Running PicardTools converter for file: ' + i + "\n")
                sample_id, sep = i.split(".bam")
                #print(sample_id)
                #picard.jar SamToFastq I=11354-L994H.bam FASTQ=11354-L994.fastq.gz
                print("picard.jar SamToFastq I=" + i + " FASTQ=" + sample_id + ".fastq.gz")
                os.system("picard.jar SamToFastq I=" + i + " FASTQ=" + sample_id + ".fastq.gz")
                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "picard.jar SamToFastq I=" + i + " FASTQ=" + sample_id + ".fastq.gz")


def changeBamChrNotation(before122 = "chr", after122 = "", beforeXYM = "X,Y,M", afterXYM = "23,24,M", suffix = "_reheader"):
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< BATCH BAM CHR NOTATION - Replace BAM chromosome notation >")
    fout.write("\n< BATCH BAM CHR NOTATION - Replace BAM chromosome notation >\n")

    for i in os.listdir(cwd):
        if re.search('.bam$', i) != None:
            print('\t> Replacing notation using samtools reheader for file: ', i)
            fout.write('\t> Replacing notation using samtools reheader for file: ' + i + "\n")
            sample_id, sep = i.split(".bam")
            
            print("samtools view -H " + i + " | sed -e 's/SN:" + before122 + "/SN:" + after122 + "/' -e 's/SN:" + beforeXYM.split(",")[0] + "/SN:" + afterXYM.split(",")[0] + "/' -e 's/SN:"  + beforeXYM.split(",")[1] + "/SN:" + afterXYM.split(",")[1] + "/' -e 's/SN:"  + beforeXYM.split(",")[2] + "/SN:" + afterXYM.split(",")[2] + "/' | samtools reheader - " + i + " > " + sample_id + suffix + ".bam")
            os.system("samtools view -H " + i + " | sed -e 's/SN:" + before122 + "/SN:" + after122 + "/' -e 's/SN:" + beforeXYM.split(",")[0] + "/SN:" + afterXYM.split(",")[0] + "/' -e 's/SN:"  + beforeXYM.split(",")[1] + "/SN:" + afterXYM.split(",")[1] + "/' -e 's/SN:"  + beforeXYM.split(",")[2] + "/SN:" + afterXYM.split(",")[2] + "/' | samtools reheader - " + i + " > " + sample_id + suffix + ".bam")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools view -H " + i + " | sed -e 's/SN:" + before122 + "/SN:" + after122 + "/' -e 's/SN:" + beforeXYM.split(",")[0] + "/SN:" + afterXYM.split(",")[0] + "/' -e 's/SN:"  + beforeXYM.split(",")[1] + "/SN:" + afterXYM.split(",")[1] + "/' -e 's/SN:"  + beforeXYM.split(",")[2] + "/SN:" + afterXYM.split(",")[2] + "/' | samtools reheader - " + i + " > " + sample_id + suffix + ".bam")

            print("samtools index " + sample_id + suffix + ".bam")
            os.system("samtools index " + sample_id + suffix + ".bam")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools index " + sample_id + suffix + ".bam")

#sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' 
#" | sed -e 's/SN:" + before122 + "/SN:" + after122 + "/' -e 's/SN:" + beforeXYM.split()[0] + "/SN:" + afterXYM.split()[0] + "/' -e 's/SN:"  + beforeXYM.split()[1] + "/SN:" + afterXYM.split()[1] + "/' -e 's/SN:"  + beforeXYM.split()[2] + "/SN:" + afterXYM.split()[2] + "/' | samtools reheader - " + i + " > " + sample_id + suffix + ".bam"


def pileupCaller(referenceGenome, posListBed, bamExtension = ".q30.rmdup.bam", minBQ = 30, minMQ = 30):
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< PILEUPER - Generate pileup from BAM files >")
    fout.write("\n< PILEUPER - Generate pileup from BAM files >\n")

    for i in os.listdir(cwd):
        if re.search((bamExtension+'$'), i) != None:
            print('\t> Gathering pileup information for: ', i)
            fout.write('\t> Gathering pileup information for: ' + i + "\n")
            sample_id, sep = i.split(bamExtension)
            print("samtools mpileup -Q " + str(minBQ) + " -q " + str(minMQ) + " -B -f " + referenceGenome + " -l " + posListBed + " " + i + " > " + sample_id + ".pileupsamtools.txt")
            os.system("samtools mpileup -Q " + str(minBQ) + " -q " + str(minMQ) + " -B -f " + referenceGenome + " -l " + posListBed + " " + i + " > " + sample_id + ".pileupsamtools.txt")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools mpileup -Q " + str(minBQ) + " -q " + str(minMQ) + " -B -f " + referenceGenome + " -l " + posListBed + " " + i + " > " + sample_id + ".pileupsamtools.txt")

def pileup2ped(mapList, toBed = True, exclTerminals = False):
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< PILEUP2PED - Generate pseudo-haploid PLINK files from pileup >")
    fout.write("\n< PILEUP2PED - Generate pseudo-haploid PLINK files from pileup >\n")

    for i in os.listdir(cwd):
        if '.pileupsamtools.txt' == i[-19:]:
            if os.path.getsize(i) > 0:
                sample_id, rest = re.split('.pileupsamtools.txt', i)

                print("awk '{print $1" + '"-"' + "$2}' " + i + " > temporary1")
                os.system("awk '{print $1" + '"-"' + "$2}' " + i + " > temporary1")
                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "awk '{print $1" + '"-"' + "$2}' " + i + " > temporary1")

                print("awk '{print $1" + '"-"' + "$4}' " + mapList + " | grep -n -x -w -F -f temporary1 | awk '{split($0, array, " + '":"' + "); print(array[1]) }' | awk 'FNR==NR{a[$1];next}FNR in a' - " + mapList + " > " + sample_id + ".map")
                os.system("awk '{print $1" + '"-"' + "$4}' " + mapList + " | grep -n -x -w -F -f temporary1 | awk '{split($0, array, " + '":"' + "); print(array[1]) }' | awk 'FNR==NR{a[$1];next}FNR in a' - " + mapList + " > " + sample_id + ".map")
                fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "awk '{print $1" + '"-"' + "$4}' " + mapList + " | grep -n -x -w -F -f temporary1 | awk '{split($0, array, " + '":"' + "); print(array[1]) }' | awk 'FNR==NR{a[$1];next}FNR in a' - " + mapList + " > " + sample_id + ".map")
                os.remove("temporary1")

                pedout = open(sample_id + '.ped', 'w')
                pileupin = open(i, 'r')     
                snpcount, discar = re.split(' ', os.popen('wc -l ' + i).read())
                snpcount = int(snpcount)
                pedout.write(sample_id + ' ' + sample_id + ' 0 0 0 1 ')
                pileupin = open(i, 'r')
                counter = 1
                countbads = 0
                for line in pileupin:
                    if counter <= snpcount:
                        line = line.rstrip()
                        if line.split("\t")[3] == "0":
                                pedout.write('0 0 ')
                                countbads += 1
                        else:
                            chro, coord, refb, cov, rbase, qual = line.split('\t')
                            counter += 1
                            rbase = rbase.replace(".",refb.capitalize())
                            rbase = rbase.replace(",",refb.capitalize())
                            if int(cov) == 1:
                                if r"*" in rbase or r"+" in rbase or r"-" in rbase:
                                    pedout.write('0 0 ')
                                elif r"$" in rbase:
                                    rbase = rbase[0]
                                    if exclTerminals == True: 
                                        pedout.write('0 0 ')
                                    else:
                                        pedout.write(rbase.capitalize() + ' ' + rbase.capitalize() + ' ')

                                elif r"^" in rbase:
                                    rbase = rbase[2]
                                    if exclTerminals == True:
                                        pedout.write('0 0 ')
                                    else:
                                        pedout.write(rbase.capitalize() + ' ' + rbase.capitalize() + ' ')
                                else:
                                    pedout.write(rbase.capitalize() + ' ' + rbase.capitalize() + ' ')
                            
                            ## If there is more than one call for this position
                            elif int(cov) > 1:
                                listBases = list(rbase)
                                ## Remove indels
                                while r"-" in listBases:
                                    for i in listBases:
                                        if i == r"-":
                                            indexi = listBases.index(r"-")
                                            delL = listBases[indexi+1]
                                            del listBases[indexi:(indexi+int(delL)+2)] # delete indel annotation
                                            del listBases[indexi-1] # delete indel
                                while r"+" in listBases:
                                    for i in listBases:
                                        if i == r"+":
                                            indexi = listBases.index(r"+")
                                            delL = listBases[indexi+1]
                                            del listBases[indexi:(indexi+int(delL)+2)] # delete indel annotation
                                            del listBases[indexi-1] # delete indel

                                ## Remove start and end of read calls is requested (^ and $)
                                while r"^" in listBases:
                                    for i in listBases:
                                        if i == r"^":
                                            indexi = listBases.index(r"^")
                                            del listBases[indexi] # delete ^
                                            del listBases[indexi] # delete the base quality after the ^
                                            if exclTerminals == True:
                                                del listBases[indexi] # delete the actual variant 
                                while r"$" in listBases:
                                    for i in listBases:
                                        if i == r"$":
                                            indexi = listBases.index(r"$")
                                            del listBases[indexi] # delete $
                                            if exclTerminals == True:
                                                del listBases[indexi-1] # delete the actual variant 

                                if len(listBases) != 0:
                                    chosenbase = random.choice(listBases)
                                    pedout.write(chosenbase.capitalize() + ' ' + chosenbase.capitalize() + ' ')
                                else:
                                    countbads += 1
                                    pedout.write('0 0 ')
                    else:
                        pass
                pedout.close()

                if toBed == True:
                    print("plink --file " + sample_id + " --make-bed --out " + sample_id + " --allow-no-sex --keep-allele-order")
                    os.system("plink --file " + sample_id + " --make-bed --out " + sample_id + " --allow-no-sex --keep-allele-order")
                    fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "plink --file " + sample_id + " --make-bed --out " + sample_id + " --allow-no-sex --keep-allele-order")
                    os.remove(sample_id + ".ped")
                    os.remove(sample_id + ".map")


def schmutzi(referenceMTgenome, schmutziFolder, threads = 4, bamExtension = ".q30.rmdup.bam"):
    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< SCHMUTZI - Estimate contamination using schmutzi >")
    fout.write("\n< SCHMUTZI - Estimate contamination using schmutzi >\n")

    for i in os.listdir(cwd):
        if re.search(bamExtension + "$", i) != None:
            print('\t> Estimating contamination using schmutzi for file: ', i)
            fout.write('\t> Estimating contamination using schmutzi for file: ' + i + "\n")
            sample_id, sep = i.split(bamExtension)

            print("samtools calmd -b " + i + " " + referenceMTgenome + " > " + sample_id + "_MD" + bamExtension)
            os.system("samtools calmd -b " + i + " " + referenceMTgenome + " > " + sample_id + "_MD" + bamExtension)
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools calmd -b " + i + " " + referenceMTgenome + " > " + sample_id + "_MD" + bamExtension)
                
            print("samtools index " + sample_id + "_MD" + bamExtension)
            os.system("samtools index " + sample_id + "_MD" + bamExtension)
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools index " + sample_id + "_MD" + bamExtension)
                
            print("contDeam.pl --uselength --library double --out " + sample_id + "_MD " + referenceMTgenome + " " + sample_id + "_MD" + bamExtension)
            os.system("contDeam.pl --uselength --library double --out " + sample_id + "_MD " + referenceMTgenome + " " +sample_id + "_MD" + bamExtension)
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "contDeam.pl --uselength --library double --out " + sample_id + "_MD" + referenceMTgenome + " " + sample_id + "_MD" + bamExtension)

            print("schmutzi.pl --notusepredC --t " + str(threads) + " --ref " + referenceMTgenome + " --out " + sample_id + "_MD_npred " + sample_id + "_MD" + schmutziFolder + "/alleleFreqMT/eurasian/freqs/ " + sample_id + "_MD" + bamExtension)
            os.system("schmutzi.pl --notusepredC --t " + str(threads) + " --ref " + referenceMTgenome + " --out " + sample_id + "_MD_npred " + sample_id + "_MD" + schmutziFolder + "/alleleFreqMT/eurasian/freqs/ " + sample_id + "_MD" + bamExtension)
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "schmutzi.pl --notusepredC --t " + str(threads) + " --ref " + referenceMTgenome + " --out " + sample_id + "_MD_npred " + sample_id + "_MD" + schmutziFolder + "/alleleFreqMT/eurasian/freqs/ " + sample_id + "_MD" + bamExtension)

            print("schmutzi.pl --t " + str(threads) + " --ref " + referenceMTgenome + " --out " + sample_id + "_MD_npred " + sample_id + "_MD " + schmutziFolder + "/alleleFreqMT/eurasian/freqs/ " + sample_id + "_MD" + bamExtension)
            os.system("schmutzi.pl --t " + str(threads) + " --ref " + referenceMTgenome + " --out " + sample_id + "_MD_npred " + sample_id + "_MD " + schmutziFolder + "/alleleFreqMT/eurasian/freqs/ " + sample_id + "_MD" + bamExtension)
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "schmutzi.pl --t " + str(threads) + " --ref " + referenceMTgenome + " --out " + sample_id + "_MD_npred " + sample_id + "_MD " + schmutziFolder + "/alleleFreqMT/eurasian/freqs/ " + sample_id + "_MD" + bamExtension)

def hapConX(hapROHfolder = "~/Software/hapROH/", minMapQ = 30, minBaseQ = 30, bamExtension = ".q30.rmdup_reheader.bam"):
    os.chdir(originalWD)
    from os import path
    cwd = os.getcwd()
    print("\n< HAPCONX - Estimate contamination using hapConX >")
    fout.write("\n< HAPCONX - Estimate contamination using hapConX >\n")

    summaryOut = open('HapConX_Results.txt', 'w')
    summaryOut.write("Sample\tContamination\t95pc_CI\tSNPs_Used\tCompiled\n")
    for i in os.listdir(cwd):
        if re.search(bamExtension + "$", i) != None:
            print('\t> Estimating contamination using hapConX for file: ', i)
            fout.write('\t> Estimating contamination using hapConX for file: ' + i + "\n")
            sample_id, sep = i.split(bamExtension)

            print("samtools mpileup --positions " + hapROHfolder + "hapCon_supportFiles/chrX_1000G.bed" + " -r X -q " + str(minMapQ) + " -Q " + str(minBaseQ) + " -o " + sample_id + ".mpileup " + i)
            os.system("samtools mpileup --positions " + hapROHfolder + "hapCon_supportFiles/chrX_1000G.bed" + " -r X -q " + str(minMapQ) + " -Q " + str(minBaseQ) + " -o " + sample_id + ".mpileup " + i)
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "samtools mpileup --positions " + hapROHfolder + "hapCon_supportFiles/chrX_1000G.bed" + " -r X -q " + str(minMapQ) + " -Q " + str(minBaseQ) + " -o " + sample_id + ".mpileup " + i)
                
            print("python3 " + hapROHfolder + "/bam/hapCONX.py -m " + sample_id + ".mpileup -r " + hapROHfolder + "hapCon_supportFiles/chrX_1000G.hdf5 --meta " + hapROHfolder + "hapCon_supportFiles/meta_df_all.csv")
            os.system("python3 " + hapROHfolder + "/bam/hapCONX.py -m " + sample_id + ".mpileup -r " + hapROHfolder + "hapCon_supportFiles/chrX_1000G.hdf5 --meta " + hapROHfolder + "hapCon_supportFiles/meta_df_all.csv")
            fout.write(time.strftime("\n\t%H:%M:%S %Z ") + "\t" + "python3 " + hapROHfolder + "/bam/hapCONX.py -m " + sample_id + ".mpileup -r " + hapROHfolder + "hapCon_supportFiles/chrX_1000G.hdf5 --meta " + hapROHfolder + "hapCon_supportFiles/meta_df_all.csv")

            if path.exists(sample_id + ".hapCon.txt") == True:
                hapConFileIn = open(sample_id + ".hapCon.txt", 'r')
                snpsHapCon = hapConFileIn.readline().split(" ")[10].rstrip()
                hapConFileIn.readline()
                hapConFileIn.readline()
                lastL = hapConFileIn.readline()
                contamHapCon, confId = lastL.split(" ")[5].rstrip(), "(" + lastL.split("(")[1].rstrip()
                compi = contamHapCon + " " + confId + " - " + snpsHapCon
                #print(snpsHapCon)
                #print(contamHapCon)
                #print(confId)
                #print(compi)
                summaryOut.write(sample_id + "\t" + contamHapCon + "\t" + confId + "\t" + snpsHapCon + "\t" + compi + "\n")
                
    summaryOut.close()


            #for line in pileupin:
            #    if counter <= snpcount:
            #        line = line.rstrip()
            #        if line.split("\t")[3] == "0":


def bamDemux(inBam,sampleList):
    """
    BAM demultiplexing. P7 is barcode1 (bc), P5 is barcode2 (b2) reverse-complemented.

    Input:
    Foo.bam

    Output:
    Sample1.fastq
    Sample2.fastq

    :param inBam:
    (str)   String with input BAM file.

    :param sampleList:
    (str)   Tab-spaced list of samples and indices in the following order: Sample\tP7\tP5reverse-complemented
    """

    os.chdir(originalWD)
    cwd = os.getcwd()
    print("\n< DEMUXQ - Demultiplexing BAM into individual sample files >")
    fout.write("\n< DEMUXQ - Demultiplexing BAM into individual sample files >")

    sampIn = open(sampleList, 'r')
    sampDict = {}

    

#    if(len(sampIn.readline().split("\t"))) == 2:
#        sampIn.seek(0)
#
#        for line in sampIn:
#            if str.split(line)[0] in sampDict.values():
#                print("\t> ERROR: Found duplicated index ",str.split(line)[0])
#                print("\t> Please ensure all index sequences are unique. Quitting now...")
#                quit()
#
#            if str.split(line)[1] in sampDict.keys():
#                print("\t> ERROR: Found duplicated sample ",str.split(line)[1])
#                print("\t> Please ensure all sample names are unique. Quitting now...")
#                quit()
#            else:
#                sampDict[str.split(line)[1]] = str.split(line)[0]
    print(len(sampIn.readline().split("\t")))
    print(sampDict)

    if len(sampIn.readline().split("\t")) == 3:
        sampIn.seek(0)

        print("AAAAA")

        for line in sampIn:
            if str.split(line)[0] in sampDict.values():
                print("\t> ERROR: Found duplicated index ",str.split(line)[0])
                print("\t> Please ensure all index sequences are unique. Quitting now...")
                quit()

            if str.split(line)[1] in sampDict.keys():
                print("\t> ERROR: Found duplicated sample ",str.split(line)[1])
                print("\t> Please ensure all sample names are unique. Quitting now...")
                quit()
            else:
                sampDict[str.split(line)[1]+"\\\+"+str.split(line)[2]] = str.split(line)[0]

    sampIn.close()
    sampList=list(sampDict.keys())

    ## Check for previous files and delete them
    #lnum = inBam.split("L00")[1][0]
    #for i in sampList:
    #    for items in os.listdir(cwd):
    #        if sampDict[i] in items:
    #            try:
    #                os.remove(items)
    #            except:
    #                pass
    #try:
    #    os.remove("Unindexed_S0_" + "L00" + lnum + "_R1_001.fastq")
    #except:
    #    pass

    ## Create results-related dicts
    sampNum = {}
    seqa = 1
    for i in sampDict.keys():
        sampNum[i] = seqa
        seqa += 1

    sampFinds = {}
    for i in sampDict.keys():
        sampFinds[i] = 0
    countUndet = 0

    ## Iterate over the BAM file
    print(sampDict)
    
    for onekey in sampDict.keys():
        print(onekey)
        print(onekey.split("\\\\+")[0])
        print(sampDict[onekey])
        
        #print(aaah)

        print("/home/dfernandes/Software/samtools-1.15/installfolder/bin/samtools view -@ 6 -d BC:" + onekey.split("\\\\+")[0] + "ATCTCG " + inBam + " -o " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_P7.bam")
        os.system("/home/dfernandes/Software/samtools-1.15/installfolder/bin/samtools view -@ 6 -d BC:" + onekey.split("\\\\+")[0] + "ATCTCG " + inBam + " -o " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_P7.bam")

        print("/home/dfernandes/Software/samtools-1.15/installfolder/bin/samtools view -@ 6 -d B2:" + onekey.split("\\\\+")[1] + "GTGTAG " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_P7.bam" + " -o " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_L001_R1_001.bam")
        os.system("/home/dfernandes/Software/samtools-1.15/installfolder/bin/samtools view -@ 6 -d B2:" + onekey.split("\\\\+")[1] + "GTGTAG " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_P7.bam" + " -o " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_L001_R1_001.bam")
        
        print("bedtools bamtofastq -i " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_L001_R1_001.bam" + " -fq " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_L001_R1_001.fastq")
        os.system("bedtools bamtofastq -i " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_L001_R1_001.bam" + " -fq " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_L001_R1_001.fastq")

        print("nohup gzip " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_L001_R1_001.fastq")
        os.system("nohup gzip " + str(sampDict[onekey]) + "_S" + str(sampNum[onekey]) + "_L001_R1_001.fastq")


def main():
    output('Output', 'NextSeq', "singleEnd")  ##Must run. ALWAYS

    ## SELECT NOW THE DESIRED FUNCTIONS TO RUN, BY (UN)COMMENTING THEM

        ## Reference Genome Handling
    # chrConcat('/home/dfernandes/NGSdata/Indexes/chroms_hg19', overwrite=True, mtChrAtEnd=True)#, deleteFaFiles=True) ##Path must be inside quotation marks. Please check 'help(chrConcat)' for more information
    # makeIndexGenome('/home/dfernandes/NGSdata/Script_projects/small_edits/Cavelana500.fasta', bwa=True, samtools=True, gatkpicard=True, overwrite=True)  ##Path must be inside quotation marks. Please check 'help(makeIndexGenome)' for more information.

        ## FASTQ Handling
    # fastqdemux("Mar22aaaaa_S0_L001_R1_001.fastq", "barcodes_list_revComp5", mismatch=0)
    #bamDemux("HFFF5DRX2_2_20220630B_20220701.bam", "indeces_missingSample_199641-2")
    # concatenateFQGZ()

        ## Sequencing Data Alignment Pipeline
    #fastqcAnalysis()
    #cutadapt(nThreads=8)
    #trimFastq(3)
    #bwaAlign(ref_gen='/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa',nThreads=8, overwrite=True)
    #tagAndConvertSb(ref_gen='/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa',group_tag="Osor-SPetar",qualThreshold=30)
    #samToolsSort(version="1.3.1")
    #fastqcBeforeRmdup()
    #rmvDups(picard=False)
    #makeIndexBam()
    #mapdamage('/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa',stats = False)
    # rescaleMapDamage('/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa')
    #fastqcAfterRmdup()

        ## Sexing Assessment
    #idxstats()
    sexing()
    #sexingEJ(single_out=True, st_out=True)

        ## Haplogroup Assessment
    #mtDNAphymer('/home/dfernandes/Software/phy-mer-master/PhyloTree_b16_k12.txt')
    # Yfitter('/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa', '/home/dfernandes/Software/YFitter_v0.2/karafet_sites_b37_tab','/home/dfernandes/Software/YFitter_v0.2/karafet_tree_b37',mergeLanes = True)
    # yleaf(refPath='/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa',refVersion='hg19')

        ## Contamination Assessment
    #angsd_contamination(".q30.rmdup.bam","R")
    # calico(refLoc="/home/dfernandes/Software/calico_0.2/human.0.95.consensus.fasta",calLoc="/home/dfernandes/Software/calico_0.2/calico.0.2.py", qualThreshold=30,multipleLanes=False)  # Possibly fixed?
    #schmutzi(referenceMTgenome = "/home/dfernandes/NGSdata/Indexes/rCRS/rCRS.fasta", schmutziFolder = "/home/dfernandes/Software/schmutzi/", threads = 8, bamExtension = ".q30.rmdup.bam")
    #changeBamChrNotation(before122 = "chr", after122 = "", beforeXYM = "X,Y,M", afterXYM = "X,Y,MT",suffix = "_reheader")
    #hapConX() #Requires reheader BAMs without "chr" in notation

        ## Summary Tests
    #coverage(contigFileLocation="/home/dfernandes/NGSdata/Indexes/chroms_hg19/contigs.txt")
    #readlengths()

        ## End Analysis (Compile Data)
    #end_analysis()
    # mergedEnd_Analysis()

        ## BAM Editing
    # reorderBAM("/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa",".q30.rmdup.bam",True)
    # editHeaderBAMbatch(extensionBAM=".q30.rmdup.bam",RGid="Houtaomuga",rename=True)
    #changeBamChrNotation(before122 = "chr", after122 = "", beforeXYM = "X,Y,M", afterXYM = "23,24,M",suffix = "_reheader")

        ## BAM Handling and Analysis
    # subsampling(0.25,5, forSexing=False)
    # extract_reads_pos(6, 396320, 396323)
    # extract_reads_seq("CTCGCTGAACCGGATCGATGTGTACTGAGATCCCCTATCCGTATGGTGGATAACGTCTTTCAGGTCGAGTACGCCTTCTTGTTGGC", 0.4, 1, oneFile="ARA1%_R1_001.trimmed.fastq")
    # hirisplex_lct("/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa")
    # mergeBamsBatch(qualThreshold=30, remSourceFiles = False)
    # mergeBams(sample_id="Thingy", qualThreshold=25, version="1.3.1")  ##Sample_id AND files to merge cannot have the word "merged" on them, because that breaks the script
    #pca_lc('/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa', '/home/dfernandes/NGSdata/Datasets/laz16/lowResPCA', '/home/dfernandes/NGSdata/Datasets/laz16/lowResPCA/dataPP.bim', '/home/dfernandes/NGSdata/Datasets/laz16/lowResPCA/laz16PPmap.list', indexGen=False, makepileup=True, minNumReadsPileup = 15000, makeped=True, pedfamily="SP", makemap=True, conv2bin=True, mergeIt=False, reduceDSsnp=False, genoPlink=0.1, reduceDSind=False, popList="/home/dfernandes/NGSdata/Datasets/laz16/lowResPCA/EU_Modern_Pops", PCAparamFile=True, plot="/home/dfernandes/NGSdata/Datasets/laz16/Populations_Laz16_IE_GR_PL_orig.csv")#, useJava8=True)
    #pca_lc('/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa', '/home/dfernandes/NGSdata/Datasets/genocalling/1240k_genocalling', '/home/dfernandes/NGSdata/Datasets/genocalling/1240k_genocalling/1240K.bim', '/home/dfernandes/NGSdata/Datasets/genocalling/1240k_genocalling/1240K_chrNums_chrXY.bim.list', indexGen=False, makepileup=True, minNumReadsPileup = 15000, makeped=True, pedfamily="Romito", makemap=True, conv2bin=True, mergeIt=False, reduceDSsnp=False, genoPlink=0.1, reduceDSind=False, popList="/home/dfernandes/NGSdata/Datasets/laz16/lowResPCA/EU_Modern_Pops", PCAparamFile=True, plot="/home/dfernandes/NGSdata/Datasets/laz16/Populations_Laz16_IE_GR_PL_orig.csv")#, useJava8=True)
    # admix(2,12,reorder=FALSE)
    # admix2(2,13, multithread_cores=8,input_file="Laz16ALLPolarizeEIG1200KAutosomesOrdered")
    # admix3(k_num=3,reps=10,multithread_cores=8,input_file="Laz16ALLPolarizeEIG1200KAutosomesOrdered")
    # imputeWithBeagle(curateRefPanel=True,bcftoolsNthreads=1)
    # batchBAMtoFASTQ(picard = True)
    #pileupCaller(bamExtension = ".q30.rmdup_reheader.bam", minBQ = 30, minMQ = 30, referenceGenome = '/home/dfernandes/NGSdata/Indexes/chroms_hg19_nochr/full_karyo.fa', posListBed = "/home/dfernandes/NGSdata/Indexes/200kwien/200kwienNoChr.bed") #"/home/dfernandes/NGSdata/Indexes/chroms_hg19/hirisplexs.bed"
    #pileup2ped(mapList = "/home/dfernandes/NGSdata/Indexes/200kwien/200K_noChr.map", exclTerminals = False)
    # hirisplexs("/home/dfernandes/NGSdata/Indexes/chroms_hg19_nochr/full_karyo.fa")

        ## Varied Experiments
    # damageRestrict(samIn="C3.q30.rmdup.bam",refGenome="/home/dfernandes/NGSdata/Indexes/chroms_hg19/full_karyo.fa")
    # CpGdamage()
    # transFst(bimFile="Yoruba.bim")
    # blastnQuery(queryLink="ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/")
    # moveFilesBatch(byExtension = "q30.rmdup.bam/q30.rmdup.bam.bai", bySample = "AC/EC/GNT/SCH", copy = True, singleFolder=False)

    exit()


if __name__ == "__main__":
    main()

# 12, 6643657, 6647536 - should be unmethylated/damaged
# 6, 25726291, 25726790 - should be methylated/damaged
# 6, 25727137, 25727573  - should be methylated/damaged

# "q30\.rmdup\.bam/q30\.rmdup\.bam\.bai"
# "AC/EC/GNT/SCH"

# TODO:
# PRIORITY: For single-threaded steps, add argument asking how many parallel commands to run. Will significantly improve speed
# Try to make a quicker mapdamage graph with my own script, based on CpGdamage/damageRestrict
# Problems with directory changing in end_analysis because of mtDNA or pca_lc. Fix with "originalWD" as global variable
# Add readlength averages and stdv in end_analysis
# Put error check for sexingEJ as it gives no error if there is no bar.chr.txt
# Check on output() if picard (or any other JAVA software can be called directly or need "java -jar location/file.jar"
# bwaAlign does not recognize seqtk files if trimBam creates them right before it. Only recognizes if seqtk files already exist when "output()" runs
# No command sor any output for trimBam



# NEXTSEQ problems:
# Summary stats do not work
# sample_ids are not correctly called in end_analysis function (probably source of all problems)
