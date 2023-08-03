#%%

import time, os, re, pandas as pd, subprocess, tempfile
from datetime import datetime
import functions as fn, searchTools as st
import configSettings as cfg
from json import loads, dumps

os.chdir(os.path.dirname(__file__))

freyja = "/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/data/ncov-ww_upto_230728.tsv"
freyja1 = fn.formatFreyjaLineage(freyja).reset_index(['resid','coverage'])
freyja2 = fn.formatFreyjaLineage(freyja,summarized=True).reset_index(['resid','coverage'])

print(freyja1)
print(freyja2)

freyja = freyja1.merge(freyja2,how="outer",on="file")

print(freyja)

freyja.to_csv("../results/ncov-ww_upto_230728.csv")






# %%
