from functions import *

os.chdir(os.path.dirname(__file__))
freyjaPath = "/nfs/APL_Genomics/virus_covid19/ww_covid/university_seq/freya_test\results"
refPath = "/nfs/APL_Genomics/virus_covid19/ww_covid/university_seq/nCoV-2019.reference"
BAMdata = getBAMdb(freyjaPath, BAMdb="./data/BAMdb.txt", BAMdata = "./data/BAMdb.csv")
freyja = runFrejya(BAMdata[BAMPathCol].values.tolist(), outDir = "./results", ref = "./data/ncov.fasta")
print(freyja)
print(formatFreyjaOutput(freyja))