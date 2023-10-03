#%%

import time, os, re, pandas as pd, subprocess, tempfile, pathlib
from datetime import datetime
import functions as fn, searchTools as st
import configSettings as cfg
from json import loads, dumps

os.chdir(os.path.dirname(__file__))

rand = "/nfs/APL_Genomics/scratch/230927_N_I_008_WW_rand_out/230927_N_I_008_WW_rand/results/freyja/demix/230927_N_I_008_WW_rand_freyja_agg.tsv"

fn.formatFreyjaLineage(rand).to_csv("/nfs/APL_Genomics/scratch/230927_N_I_008_WW_rand_out/230927_N_I_008_WW_rand/results/freyja/demix/230927_N_I_008_WW_rand_freyja_agg_out.csv")

rand = "/nfs/APL_Genomics/scratch/230927_N_I_008_WW_seq_out/230927_N_I_008_WW_seq/results/freyja/demix/230927_N_I_008_WW_seq_freyja_agg.tsv"

fn.formatFreyjaLineage(rand).to_csv("/nfs/APL_Genomics/scratch/230927_N_I_008_WW_seq_out/230927_N_I_008_WW_seq/results/freyja/demix/230927_N_I_008_WW_seq_freyja_agg_out.csv")

# tsv = "/nfs/APL_Genomics/virus_covid19/ww_covid/cumulative output/230615_N_I_059_out.tsv"
# csv = "/nfs/APL_Genomics/virus_covid19/ww_covid/cumulative output/230615_N_I_059_out.csv"
# xlsx = "/nfs/APL_Genomics/virus_covid19/ww_covid/cumulative output/230615_N_I_059_out.xlsx"

# dir = "/nfs/APL_Genomics/virus_covid19/ww_covid/cumulative output/data/"
# ids = pd.read_excel("/nfs/APL_Genomics/virus_covid19/ww_covid/cumulative output/cumulative_ww_out.xlsx")

# def absoluteFilePaths(directory):
#     for dirpath,_,filenames in os.walk(directory):
#         for f in filenames:
#             yield os.path.abspath(os.path.join(dirpath, f))

# files = absoluteFilePaths(dir)

# dfs = [st.importToDataFrame(file) for file in files]
# df_siteID = [df for df in dfs if "SiteID" in df.columns]
# df_siteID = pd.concat(df_siteID)

# df_accession = [df for df in dfs if "Sample" in df.columns]
# df_accession = pd.concat(df_accession)

# ids = ids.merge(df_siteID, on="SiteID", how = "outer")
# ids = ids.merge(df_accession, on="Sample", how = "outer")

# print(ids)

# ids.to_excel("/nfs/APL_Genomics/virus_covid19/ww_covid/cumulative output/cumulative_ww_out_all.xlsx")













# %%
