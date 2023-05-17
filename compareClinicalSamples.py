from functions import *

os.chdir(os.path.dirname(__file__))
seqPath = "/nfs/APL_Genomics/virus_covid19/routine_seq/"
metadata = "/nfs/Genomics_DEV/projects/nextstrain/EBS-parser/results/allData.csv"
BAMdb = "/nfs/Genomics_DEV/projects/nextstrain/EBS-parser/data/routineSeqDB.txt"
BAMfiles = "./data/BAMfiles.csv"
variantsOut = "./results/variants.tsv" 
depthsOut = "./results/depth.out"
ref = "./data/ncov.fasta"
workDir = "./work"
synBAMs = "./results/synBAMs.bam"
synPropOut = "./results/synProp.tsv"
frejyaOut = "./results/freyja.tsv"
compareOut = "./results/results.tsv"

BAMs = getBAMdb(seqPath = seqPath, metadata = metadata, BAMdb = BAMdb, BAMfiles = BAMfiles)