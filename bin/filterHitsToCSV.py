#!/usr/bin/env python
import argparse
import json, gzip, sys
import pandas as pd
def main():
    parser = argparse.ArgumentParser(description='gets results field from jason output file of Famli')
    parser.add_argument('--jasonH', '-j', help='input jason file of hits', required=True)
    parser.add_argument('--outF', '-o', help='name of output CSV file', default="output.fasta")
    args = parser.parse_args()
    hitsInJason = args.jasonH
    outputF = args.outF
    try:
        fh = gzip.open(hitsInJason,"rt")
        sampleJson = json.load(fh)
        df = pd.DataFrame(sampleJson['results'])
    except: 
        fh = open(hitsInJason,"rt")
        sampleJson = json.load(fh)
        df = pd.DataFrame(sampleJson)

    df.to_csv(outputF)

if __name__ == '__main__':
    main()
