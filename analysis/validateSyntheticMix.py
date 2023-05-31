from functions import *

os.chdir(os.path.dirname(__file__))
seqPath = "/nfs/APL_Genomics/virus_covid19/routine_seq/2023_01_Runs"
metadata = "/nfs/Genomics_DEV/projects/nextstrain/EBS-parser/results/allData.csv"
BAMdb = "./data/BAMdb.txt"
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
synList = getSyntheticList(BAMs, n=10)
synProp = getSyntheticLineageProportions(synList)
synReads = generateSyntheticReads(BAMs = synList, outFile = synBAMs, depth = 10)
frejya = runFrejya(BAMfile = synBAMs, variantOut = variantsOut, depthsOut = depthsOut, ref = ref, outFile = frejyaOut)

#synProp = pd.read_csv(synPropOut, sep="\t")  
results = compareFreyja(frejyaOut, synProp)
print(results)
results.to_csv(compareOut, sep="\t", index = False)