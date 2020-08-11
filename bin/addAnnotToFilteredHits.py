#!/usr/bin/env python
import re
import argparse
import pandas as pd
from Bio import SeqIO


p_taxid=re.compile(r"""
(TaxID)   #the keys group
(=)  #the values starter
([^ ]+) #the values
( )  #the values terminator
""", re.VERBOSE)

p_repid=re.compile("""
        (RepID)
        (=)
        (.*)
        ($)
""", re.VERBOSE)

p_enzyme=re.compile("""
( )
(.*)
(n=)
""",re.VERBOSE)


def main():
    parser = argparse.ArgumentParser(description='gets results field from jason output file of Famli')
    parser.add_argument('--fltrdHitsCSV', '-i', help='input CSV file of filtered hits', required=True)
    parser.add_argument('--fltrdHitsFAA', '-f', help='input FAA file of filtered hits', required=True)
    parser.add_argument('--outF', '-o', help='name of output CSV file', default="output.fasta")
    args = parser.parse_args()
    #addAnnotToFilteredHits.py -i ${hitsInCSV} -f ${hitsInFaa} -o hits_in_${sname}_annot.csv
    hitsInCSV = args.fltrdHitsCSV
    hitsInFaa = args.fltrdHitsFAA
    outputF = args.outF

    df=pd.read_csv(hitsInCSV, index_col=0)
    fh=open(hitsInFaa,"rt")
    recs=SeqIO.parse(fh,'fasta')
    df.set_index('id', inplace=True)
    for rec in recs:
        if df.index.isin([rec.id]).sum()>0:
            #df.loc[rec.id,'description']=rec.description
            df.loc[rec.id,'taxid']=p_taxid.findall(rec.description)[0][2]
            df.loc[rec.id,'gene']=p_repid.findall(rec.description)[0][2].split("_")[0].lower()
            df.loc[rec.id,'description']= " ".join(p_enzyme.findall(rec.description)[0][1].split(" ")[1:])

    df.to_csv(outputF)

if __name__ == '__main__':
    main()
