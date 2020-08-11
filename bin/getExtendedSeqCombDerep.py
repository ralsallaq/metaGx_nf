#!/usr/bin/env python
from __future__ import print_function
import os, sys, time
import argparse
from Bio import SeqIO, Entrez, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import re
import pandas as pd
import numpy as np
from Bio import Entrez
import gzip

def main():
    parser = argparse.ArgumentParser(description='combining fasta files for proteins correspnding to an enzyme from uniprot and extended search from uniref90')
    parser.add_argument('--seedF', '-s', help='input seed fasta file for the enzyme (usually reviewed sequences from Uniprot)', required=True)
    parser.add_argument('--dbAA', '-db', help='input AA database fasta file (e.g. uniref90.fasta.gz)', required=True)
    parser.add_argument('--blastT', '-bt', help='input blast-like alignment table of matched records from AA db', required=True)
    parser.add_argument('--outF', '-o', help='name of output fasta file', default="output.fasta")
    args = parser.parse_args()

    seedF = args.seedF
    dbAA = args.dbAA
    outF = args.outF
    
    #blastT=pd.read_csv(args.blastT, sep="\\t", engine='c',skipinitialspace=True, header=None) 
    blastT=pd.read_csv(args.blastT, sep="\t", engine='c',skipinitialspace=True, header=None) 
    print(blastT.head())
    blastT.columns='qseqid sseqid qlen slen qstart qend sstart send gapopen mismatch pident evalue bitscore'.split(" ")
    blastT.loc[:,'sseqid']=blastT['sseqid'].apply(lambda r:r.strip())

    with open(outF,'wt') as outh:
        with open(seedF,'rt') as seedh:
            seedTtl=[]
            seedSeqs=[]
            #for rec in SeqIO.parse(seedh,"fasta"):
            for title, sequence in SimpleFastaParser(seedh):
                seedTtl.append(title)
                seedSeqs.append(sequence)
        with gzip.open(dbAA, 'rt') as dbAAh:
            for title, sequence in SimpleFastaParser(dbAAh):
                if not title.split(None,1)[0].strip() in blastT['sseqid'].values:
                    continue
                #implicit else
                #extend seed if not redundant sequence
                if not sequence in seedSeqs:
                    seedTtl.append(title)
                    seedSeqs.append(sequence)

        
        for ttl, seq in zip(seedTtl, seedSeqs):
            outh.write(">%s\n%s\n" % (ttl, seq))


if __name__ == '__main__':
    main()

