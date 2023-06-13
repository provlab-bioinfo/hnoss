import time, os, re, pandas as pd, subprocess, tempfile
from datetime import datetime
import functions as fn

os.chdir(os.path.dirname(__file__))

#file = "/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/data/GP-230403_S15_L001_R2_001.tsv"
file = "../data/results_PRL_agg.tsv"

freyja = fn.formatFreyjaLineage(file)

# print(list(freyja.columns.values.tolist()))

# print(freyja)

print(fn.collapseFreyjaLineage(freyja,['B.1.1','BA.5.2.18']))

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

