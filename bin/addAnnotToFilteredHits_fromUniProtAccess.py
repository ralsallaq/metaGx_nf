#!/usr/bin/env python
import requests, sys, os
import xml.etree.ElementTree as ET
import argparse
import numpy as np
import pandas as pd
#from io import StringIO
from io import BytesIO as StringIO
import random
import string
import json
#from datetime import datetime
#setting the random seed with system time <--- might be dangerous in parallel jobs which will run at the same time
#random.seed(datetime.now())
# create a seed from a random
seedValue = random.randrange(sys.maxsize)
# Now, Seed the random number generator
random.seed(seedValue)
print("Seed was:", seedValue)

def get_random_string(length):
    # Random string with the combination of lower and upper case
    letters = string.ascii_letters
    result_str = ''.join(random.choice(letters) for i in range(length))
    #print("Random string is:", result_str)
    return  result_str

def get_annots(requestURL):
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    
    if not r.ok:
        #r.raise_for_status()
        #sys.exit()
        pass
    
    responseBody = r.text
    json_obj = json.loads(responseBody)
    
    dataFrame = pd.DataFrame(json_obj) #columns=Index(['accession', 'comments', 'dbReferences', 'features', 'gene', 'id','info', 'keywords', 'organism', 'protein', 'proteinExistence','references', 'sequence'],dtype='object')

    accessions=[]; proteins=[]; genes=[]; ecs=set(); taxids=[]; organisms=[]; lineages=[]

    for col, rowV in dataFrame.iteritems():
        if col=='accession':
            accessions.append(rowV.iloc[0])
        if col=='organism':
            taxids.append(str(rowV.iloc[0]['taxonomy']))
            organisms.append(rowV.iloc[0]['names'][0]['value'])
            lineages.append("|".join(rowV.iloc[0]['lineage']))
        if col=='protein':
            #print(rowV.iloc[0])
            if 'recommendedName' in rowV.iloc[0].keys():
                proteins.append(rowV.iloc[0]['recommendedName']['fullName']['value'])
            elif 'submittedName' in rowV.iloc[0].keys():
                ecN=set()
                for pp in rowV.iloc[0]['submittedName']:
                    proteins.append(pp['fullName']['value'])
                    if "ecNumber" in pp.keys():
                        for eN in pp['ecNumber']:
                            ecN=ecN.union([eN['value']])
                #proteins.append(prot)
                if ecN:
                    ecs=ecs.union(ecN)
        if col=='gene':
            gene={'name':np.nan,'synonyms':[],'orfNames':[]}
            for g in rowV.iloc[0]:# going through elements of a list of dictionaries for genes
                if 'name' in g.keys():
                    gene['name']=g['name']['value']
                if 'orfNames' in g.keys():
                    for orf in g['orfNames']:
                        gene['orfNames'].append(orf['value'])
                if 'synonyms' in g.keys():
                    for syn in g['synonyms']:
                        gene['synonyms'].append(syn['value'])
            genes.append(gene)
        if col=='comment':
            for cmnt in rowV.iloc[0]:
                if cmnt['type'] == 'CATALYTIC_ACTIVITY':
                    ecs=ecs.union([cmnt['reaction']['ecNumber']])
        
    
    return accessions, list(ecs), genes, proteins, taxids, organisms, lineages

def main():
    baseURL = "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&"
    fields = 'coverage,depth,id,length,nreads,std'.split(",")
    df=pd.read_csv(sys.stdin, sep=',', header=None, index_col=0)
    df.columns = fields
    df.index.name = 'index'
    tmpOutFile = "out_"+get_random_string(8)+'.csv'

    for i, row in df.iterrows(): 
        accessNo = row['id'].split("_")[1].lower()
        requestURL=baseURL+"accession="+accessNo
        accessions, ecs, genes, proteins, taxids, organisms, lineages = get_annots(requestURL) 
        print(row['id'],organisms, taxids, accessions, ecs, genes, proteins, lineages)
        if len(accessions)>0:
            assert(accessNo==accessions[0].lower()),"accessions differ!! {} !={}".format(accessNo, accessions[0].lower())
            df.loc[i,'EC'] = ";".join(ecs)
            for g in genes:
                for k in g.keys():
                    if type(g[k]) is list:
                        df.loc[i,k] = ";".join(g[k])
                    else:
                        df.loc[i,k] = g[k]

            df.loc[i,'protein'] = ";".join(proteins)
            df.loc[i,'taxid'] = ";".join(taxids)
            df.loc[i,'organism'] = ";".join(organisms)
            df.loc[i,'lineage'] = lineages[0]

    df.to_csv(tmpOutFile)

    """
    output = StringIO()
    df.to_csv(output)
    print(output.getvalue())
    """
if __name__ == '__main__':
    main()
