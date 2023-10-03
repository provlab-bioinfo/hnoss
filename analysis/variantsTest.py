#%%

import time, os, re, pandas as pd, subprocess, tempfile
from datetime import datetime
import functions as fn, searchTools as st
import configSettings as cfg
from json import loads, dumps

os.chdir(os.path.dirname(__file__))

# files = ['../data/230810_N_I_001_WW_1_S1_L001_variants.tsv','../data/230810_N_I_001_WW_2_S2_L001_variants.tsv','../data/230810_N_I_001_WW_3_S3_L001_variants.tsv']
# variants = fn.readFreyjaVariants(files)



# print(muts)



freyja = fn.formatFreyjaLineage("../data/230810_N_I_001_WW_freyja_agg_renamed.tsv", summarized=False)
# print(freyja)
# freyja.to_excel("../results/230810_N_I_001_WW_freyja_agg_renamed.xlsx")
# freyja = fn.filterFreyjaLineage(freyja)
# fn.codeMissingAsOther(freyja).to_excel("../results/230810_N_I_001_WW_freyja_agg_renamed_filtered.xlsx")
# print(fn.normalizeValues(freyja, max = 100))
# print(fn.codeMissingAsOther(freyja, target = 1))

files = ["/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_1_S1_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_2_S2_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_3_S3_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_4_S4_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_5_S5_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_6_S6_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_7_S7_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_8_S8_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_9_S9_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_10_S10_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_11_S11_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_12_S12_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_13_S13_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_14_S14_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_15_S15_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_16_S16_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_17_S17_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_18_S18_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_19_S19_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_20_S20_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_21_S21_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_22_S22_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_23_S23_L001_variants.tsv","/nfs/APL_Genomics/virus_covid19/ww_covid/230830_N_I_002_WW/results/freyja/230830_N_I_002_WW_24_S24_L001_variants.tsv"]
mutlist = pd.read_csv('/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/data/BA.2.86-muts.txt', sep="\t", index_col = None)

variants = fn.readFreyjaVariants(files)
mutations = fn.findMutations(variants, mutlist)
count = mutations.groupby('file').size()
print(count)
print(f"BA.2.86 mutations: {len(mutlist.index)}")




# mutations.to_excel("BA.2.86-found-mutations.xlsx")

# freyja2 = fn.normalizeValues(freyja)
# print(freyja2)
# print(freyja2.sum(axis=1))
# print(fn.normalizeValues(freyja))
# print(fn.codeMissingAsOther(freyja))







# %%
