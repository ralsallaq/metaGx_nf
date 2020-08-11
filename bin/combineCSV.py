#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import glob


def main():
    files = glob.glob("./out_*.csv")
    outputFile=sys.argv[1]
    dfs = [pd.read_csv(f,index_col=0) for f in files]
    df_all = pd.concat(dfs, axis=0, ignore_index=True)
    df_all.to_csv(outputFile)
    

if __name__ == '__main__':
    main()
