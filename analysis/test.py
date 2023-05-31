import time, os, re, pandas as pd, subprocess, tempfile
from datetime import datetime
from searchTools import *
from functions import *

#region: pandas static vars
fileCol = "file"
lineageCol = "lineages"
summarizedCol = "summarized"
parentLineageCol = "parent_lineages"
abundCol = "abundances"
dominantCol = "Dominant"
dominantPercentCol = "Dominant %"
sdCol = "SD"
locationCol = "Location"
collectDateCol = "Collection Date"
residualCol = "resid"
coverageCol = "coverage"
BAMPathCol = "BAMPath"
#endregion

#file = "/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/data/GP-230403_S15_L001_R2_001.tsv"
file = "./data/results_PRL_agg.tsv"

freyja = getFreyjaLineageProportions(file)
print(freyja)
