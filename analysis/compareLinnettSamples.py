#%%
import os, functions as fn

os.chdir(os.path.dirname(__file__))

# Compare our and Linnetts samples
freyja1 = fn.formatFreyjaLineage("/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/data/compare_linnett.tsv",summarized=True)
freyja2 = fn.formatFreyjaLineage("/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/data/compare_PRL.tsv",summarized=True)

fn.compareRuns(freyja1, freyja2, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="scatter", outFile = "scatter.png")
fn.compareRuns(freyja1, freyja2, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="tukey",outFile = "tukey.png")
# %%