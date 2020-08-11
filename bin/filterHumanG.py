#!/usr/bin/env python
from __future__ import print_function
import argparse
from Bio import SeqIO, Entrez, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import re
import os, sys
import time as timeM
import pandas as pd
import numpy as np
from Bio import Entrez
import gzip

def main():

    args_parser = argparse.ArgumentParser()
    args_parser.add_argument('--input_fasta','-i', help='unfiltered human genome fasta', required=True)
    args_parser.add_argument('--output_fasta','-o', help='filtered human genome fasta', required=True)
    args_parser.add_argument('--removeID','-fltr', help='what to remove, options "EBV" or "MT" or "both"', default='both')
    args = args_parser.parse_args()

    removeID = args.removeID
    def choosefunc(removeID):
        if removeID == 'both':
            #remove Epstien-Barr virus contigs and Mitochondrial contigs
            return lambda sr: None if (sr.id == 'chrEBV') or (sr.id == 'chrM') else 1
        elif removeID == 'EBV':
            return lambda sr: None if sr.id == 'chrEBV' else 1
        elif removeID == 'MT':
            return lambda sr: None if sr.id == 'chrM' else 1
        else:
            return None

    #select function based on filter
    func = choosefunc(removeID)

    remaining=[] 
    with open(args.output_fasta,"wt") as outh:
        with gzip.open(args.input_fasta,"rt") as inh:
            for sr in SeqIO.parse(inh, format='fasta'):
                print(sr.id)
                ret=func(sr)
                if ret:
                    remaining.append(sr)
        SeqIO.write(remaining, outh,'fasta')


if __name__ == "__main__":
    main()
