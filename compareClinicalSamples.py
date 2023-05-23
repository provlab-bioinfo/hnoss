from functions import *

os.chdir(os.path.dirname(__file__))
#freyjaPath = "/nfs/APL_Genomics/virus_covid19/ww_covid/university_seq/freya_test/results/all_data"
freyjaPath = "/nfs/Genomics_DEV/projects/alindsay/Projects/wwCOV/data"
#refPath = "/nfs/APL_Genomics/virus_covid19/ww_covid/university_seq/nCoV-2019.reference"
#refPath = "./data/ncov.fasta"
refPath = "./data/nCoV-2019.reference.fasta"
BAMdata = getBAMdb(freyjaPath, BAMdb="./data/BAMdb.txt", BAMdata = "./data/BAMdb.csv")
freyja = runFrejya(BAMdata[BAMPathCol].values.tolist()[0:3], outDir = "./results", ref = refPath, refname="MN908947.3")
print(formatFreyjaOutput(freyja))