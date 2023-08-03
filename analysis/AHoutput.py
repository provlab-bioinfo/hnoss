import os, functions as fn
os.chdir(os.path.dirname(__file__))

freyja = "/nfs/APL_Genomics/apps/development/ncov-wastewater/ncov-ww_upto_230728.tsv"
freyja = fn.formatFreyjaLineage(freyja)
freyja = fn.filterFreyjaLineage(freyja, cutoff = 0.05)
freyja.to_excel("freyja_out.xlsx")


