#%%

import time, os, re, pandas as pd, subprocess, tempfile
from datetime import datetime
import functions as fn, searchTools as st
import configSettings as cfg
from pathlib import Path
from json import loads, dumps

os.chdir(os.path.dirname(__file__))

filter = True

file1 = "/nfs/APL_Genomics/scratch/compare/rand_1.0.xlsx"
file2 = "/nfs/APL_Genomics/scratch/compare/rand_0.1.xlsx"
xlab = Path(file1).stem
ylab = Path(file2).stem

def formatMe(file, summarized = True, filter = False):
    freyja = fn.formatFreyjaLineage(file, summarized=summarized)
    if filter: freyja = fn.filterFreyjaLineage(freyja,0.01)
    freyja = freyja.reset_index(level='file')
    freyja['file'] = freyja['file'].replace('[0-9]*\.[0-9]+_L001', 'L001', regex=True)
    freyja = freyja.set_index('file',append=True)
    freyja = freyja*100 + 1
    return freyja

def compare(freyja1, freyja2, xlab, ylab, log = True):
    fn.compareRuns(freyja1, freyja2, xlab=xlab, ylab=ylab, type="scatter", log=log)#, outFile = "scatter.png")
    fn.compareRuns(freyja1, freyja2, xlab=xlab, ylab=ylab, type="tukey")#,outFile = "tukey.png")

# print("Overall Summarized Lineages")
# freyja1 = formatMe(file1, summarized=True)
# freyja2 = formatMe(file2, summarized=True)

# compare(freyja1, freyja2, xlab, ylab)

print("Overall Raw Lineages")
freyja1 = formatMe(file1, summarized=False)
freyja2 = formatMe(file2, summarized=False)
compare(freyja1, freyja2, xlab, ylab, log = False)

print("Overall Raw Lineages (log-log)")
freyja1 = formatMe(file1, summarized=False)
freyja2 = formatMe(file2, summarized=False)
compare(freyja1, freyja2, xlab, ylab, log = True)


freyja1.dropna(how='all', axis=1).to_excel(Path(file1).with_name(Path(file1).stem + '_out' + Path(file1).suffix))
freyja2.dropna(how='all', axis=1).to_excel(Path(file2).with_name(Path(file2).stem + '_out' + Path(file2).suffix))


# %%
