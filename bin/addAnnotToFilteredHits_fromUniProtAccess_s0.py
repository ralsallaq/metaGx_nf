#!/usr/bin/env python
import requests, sys, os
import xml.etree.ElementTree as ET
import argparse
import numpy as np
import pandas as pd

def get_annots(requestURL):
    r = requests.get(requestURL, headers={ "Accept" : "application/xml"})
    
    if not r.ok:
        #r.raise_for_status()
        #sys.exit()
        pass
    
    responseBody = r.text
    with open("temp.xml","wt") as xmlfh:
        xmlfh.write(responseBody)
    root = ET.parse("temp.xml").getroot()
    
    
    accessions=[]; proteins=[]; genes=[]; ecs=[]; taxids=[]; organisms=[]
    for entry in root.findall("{http://uniprot.org/uniprot}entry"):
        gene=dict()
        for child in entry.getchildren():
            if child.tag=='{http://uniprot.org/uniprot}accession':
                accessions.append(child.text)
            if child.tag=='{http://uniprot.org/uniprot}gene':
                for cc in  child.getchildren():
                    gene[cc.attrib['type']]=cc.text
                genes.append(gene)
            if child.tag=='{http://uniprot.org/uniprot}protein':
                for cc in child.getchildren():
                    if cc.tag=='{http://uniprot.org/uniprot}recommendedName':
                        for c3 in cc.getchildren():
                            if c3.tag=='{http://uniprot.org/uniprot}fullName':
                                proteins.append(c3.text)
            if child.tag=='{http://uniprot.org/uniprot}comment':
                for cc in child.getchildren():
                    if cc.tag=='{http://uniprot.org/uniprot}reaction':
                        for c3 in cc.getchildren():
                            if c3.tag=='{http://uniprot.org/uniprot}dbReference' and c3.attrib['type']=='EC':
                                ecs.append(c3.attrib['id'])
            if child.tag=='{http://uniprot.org/uniprot}organism':
                for cc in child.getchildren():
                    if cc.tag=='{http://uniprot.org/uniprot}name':
                        organisms.append(cc.text)
                    if cc.tag=='{http://uniprot.org/uniprot}dbReference':
                        taxids.append(cc.attrib['id'])
                                
    #print(organisms, taxids, accessions, ecs, genes, proteins)
    return accessions, ecs, genes, proteins, taxids, organisms

def main():
    parser = argparse.ArgumentParser(description='gets results field from jason output file of Famli')
    parser.add_argument('--fltrdHitsCSV', '-i', help='input CSV file of filtered hits', required=True)
    parser.add_argument('--outF', '-o', help='name of output CSV file', default="output.fasta")
    args = parser.parse_args()
    #addAnnotToFilteredHits.py -i ${hitsInCSV} -f ${hitsInFaa} -o hits_in_${sname}_annot.csv
    hitsInCSV = args.fltrdHitsCSV
    outputF = args.outF
    df=pd.read_csv(hitsInCSV)
    print(df.head())
    baseURL = "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&"
    for i, row in df.iterrows(): 
        accessNo = row['id'].split("_")[1].lower()
        requestURL=baseURL+"accession="+accessNo

        accessions, ecs, genes, proteins, taxids, organisms = get_annots(requestURL) 
        print(organisms, taxids, accessions, ecs, genes, proteins)
        if len(accessions)>0:
            assert(accessNo==accessions[0].lower()),"accessions differ!!"
            df.loc[i,'EC'] = ",".join(ecs)
            for g in genes:
                for k in g.keys():
                    df.loc[i,k] = g[k]
            df.loc[:,'protein'] = ",".join(proteins)
            df.loc[:,'taxid'] = ",".join(taxids)
            df.loc[:,'organism'] = ",".join(organisms)

    df.to_csv(outputF)

if __name__ == '__main__':
    main()
