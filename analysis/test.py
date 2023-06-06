import time, os, re, pandas as pd, subprocess, tempfile
from datetime import datetime
from searchTools import *
from functions import *

os.chdir(os.path.dirname(__file__))
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
file = "../data/results_PRL_agg.tsv"

freyja = formatFreyjaOutput(file)

# print(list(freyja.columns.values.tolist()))

# print(freyja)

print(collapseFreyjaLineage(freyja,['B.1.1','BA.5.2.18']))

# strain = "XBB"

# aliasor = Aliasor()
# print("strain: " + strain)
# print("parent: " + aliasor.parent(strain))
# print("uncompress: " + aliasor.uncompress(strain))
# print("compress: " + aliasor.compress(strain))
# print("partial compress: " + aliasor.partial_compress(strain,up_to=50,accepted_aliases=["B","A"]))


# print(freyja)

# print(list(freyja.columns.values.tolist()))
# freyja.to_csv("../results/test.tsv",sep="\t")

