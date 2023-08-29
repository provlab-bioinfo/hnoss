#%%

import time, os, re, pandas as pd, subprocess, tempfile
from datetime import datetime
import functions as fn, searchTools as st
import configSettings as cfg
from json import loads, dumps

os.chdir(os.path.dirname(__file__))

filter = True

print("Overall Summarized Lineages")
freyja_us = fn.formatFreyjaLineage("../data/230810_N_I_001_WW_freyja_agg_renamed.tsv", summarized=True)
freyja_lin = fn.formatFreyjaLineage("../data/LIB_230815.tsv", summarized=True)

if (filter): freyja_us = fn.filterFreyjaLineage(freyja_us)
if (filter): freyja_lin = fn.filterFreyjaLineage(freyja_lin)

fn.compareRuns(freyja_us, freyja_lin, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="scatter")#, outFile = "scatter.png")
fn.compareRuns(freyja_us, freyja_lin, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="tukey")#,outFile = "tukey.png")

print("Overall Raw Lineages")
freyja_us = fn.formatFreyjaLineage("../data/230810_N_I_001_WW_freyja_agg_renamed.tsv", summarized=False)
freyja_lin = fn.formatFreyjaLineage("../data/LIB_230815.tsv", summarized=False)

if (filter): freyja_us = fn.filterFreyjaLineage(freyja_us)
if (filter): freyja_lin = fn.filterFreyjaLineage(freyja_lin)

fn.compareRuns(freyja_us, freyja_lin, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="scatter")#, outFile = "scatter.png")
fn.compareRuns(freyja_us, freyja_lin, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="tukey")#,outFile = "tukey.png")

print("Synthetic Summarized Lineages")
freyja_us = fn.formatFreyjaLineage("../data/230810_N_I_001_WW_freyja_agg_renamed.tsv", summarized=True).head(5)
freyja_lin = fn.formatFreyjaLineage("../data/LIB_230815.tsv", summarized=True).head(5)

if (filter): freyja_us = fn.filterFreyjaLineage(freyja_us)
if (filter): freyja_lin = fn.filterFreyjaLineage(freyja_lin)

fn.compareRuns(freyja_us, freyja_lin, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="scatter")#, outFile = "scatter.png")
fn.compareRuns(freyja_us, freyja_lin, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="tukey")#,outFile = "tukey.png")

print("Synthetic Raw Lineages")
freyja_us = fn.formatFreyjaLineage("../data/230810_N_I_001_WW_freyja_agg_renamed.tsv", summarized=False).head(5)
freyja_lin = fn.formatFreyjaLineage("../data/LIB_230815.tsv", summarized=False).head(5)

if (filter): freyja_us = fn.filterFreyjaLineage(freyja_us)
if (filter): freyja_lin = fn.filterFreyjaLineage(freyja_lin)

freyja_us.dropna(how='all', axis=1).to_excel("../results/230810_freyja_us.xlsx")
freyja_lin.dropna(how='all', axis=1).to_excel("../results/230810_freyja_lin.xlsx")

fn.compareRuns(freyja_us, freyja_lin, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="scatter")#, outFile = "scatter.png")
fn.compareRuns(freyja_us, freyja_lin, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="tukey")#,outFile = "tukey.png")

# %%
