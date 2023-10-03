#%%

import time, os, re, pandas as pd, subprocess, tempfile
from datetime import datetime
import functions as fn, searchTools as st
import configSettings as cfg
from json import loads, dumps

os.chdir(os.path.dirname(__file__))

# allData = pd.read_csv("/nfs/Genomics_DEV/projects/nextstrain/EBS-parser/data/BN_covid_db_export_20230627.csv", encoding_errors="ignore")#"/nfs/Genomics_DEV/projects/nextstrain/EBS-parser/results/alldata_all.csv")
# freyja_sum = pd.read_csv("/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/results/clinical_summarized.csv", encoding_errors="ignore")
# freyja_line = pd.read_csv("/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/results/clinical_lineages.csv", encoding_errors="ignore")
# freyja = freyja_sum.merge(freyja_line,how="outer",on="file")
# freyja["Key"] = freyja["file"].transform(lambda x: x.replace("_variants.tsv","").split(".")[0])

# print(freyja)
# print(allData)

# freyja = freyja.merge(allData,how="left",on="Key")

# freyja.to_csv("../results/mixed_samples_big_BN.csv")

freyja = pd.read_csv("../results/mixed_samples_big_BN.csv")

fastas = pd.DataFrame()
fastas['fastaPath'] = st.searchFlatFileDB("/nfs/Genomics_DEV/projects/nextstrain/EBS-parser/230713_fastas.txt", includeTerms=freyja['fasta'].values.tolist(), excludeTerms=['work','preconsensus'])
fastas['Key'] = fastas['fastaPath'].transform(lambda path: os.path.basename(path).rsplit('.')[0])
weights = fastas['fastaPath'].transform(lambda path: 1000000000 if bool(re.search('consensus', path)) else 1)
fastas = fastas.groupby('Key').sample(weights = weights.tolist()).reset_index()

# fastas.to_csv("../results/mixed_samples_fastas.csv")

print(fastas)

freyja = freyja.merge(fastas,how="left",on="Key")

freyja.to_csv("../results/mixed_samples_big_BN2.csv")


# ##### Aggregate the clinical samples
# freyja = fn.formatFreyjaLineage("/nfs/APL_Genomics/virus_covid19/ww_covid/all_bam/demix/clinical_agg.tsv",summarized=True)
# freyja.to_csv("../results/clinical_summarized.csv")

# freyja = fn.formatFreyjaLineage("/nfs/APL_Genomics/virus_covid19/ww_covid/all_bam/demix/clinical_agg.tsv",summarized=False)
# freyja.to_csv("../results/clinical_lineages.csv")




# #### Compare our and Linnetts samples
# freyja1 = fn.formatFreyjaLineage("/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/data/compare_linnett.tsv",summarized=True)
# freyja2 = fn.formatFreyjaLineage("/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/data/compare_PRL.tsv",summarized=True)
# freyja = fn.locateFreyjaSamples(freyja)

# fn.compareRuns(freyja1, freyja2, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", log=False)
# fn.compareRuns(freyja1, freyja2, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", log=True)
# fn.compareRuns(freyja1, freyja2, xlab="PRL COVID WW Lineage %", ylab="UofA COVID WW Lineage %", type="tukey",log=True)







# %%
