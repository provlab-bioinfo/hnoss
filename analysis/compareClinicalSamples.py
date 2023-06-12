from analysis.functions import *

os.chdir(os.path.dirname(__file__))
freyjaPath = "/nfs/APL_Genomics/virus_covid19/ww_covid/university_seq/freyja_test/results/all_data"
#freyjaPath = "/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/data"
#refPath = "/nfs/APL_Genomics/virus_covid19/ww_covid/university_seq/nCoV-2019.reference"
#refPath = "./data/ncov.fasta"
refPath = "../data/nCoV-2019.reference.fasta"
BAMdb = "../data/BAMdb.txt"
BAMdata = "../data/BAMdata.csv"

BAMdata = getBAMdb(freyjaPath, BAMdb=BAMdb, BAMdata = BAMdata)
freyja = runFrejya(BAMdata[BAMPathCol].values.tolist(), outDir = "./results", ref = refPath, refname="MN908947.3")
freyjaAgg = aggregateFreyja("./results/","./results/freyjaAgg.tsv")
plotFreyja(freyjaAgg,"freyja.png")