import os
import functions as fn
import syntheticMix as syn

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

nMixes = 10
samplesPerMix = 10
depth = 10
out = "./results/testsyn"

def createMix(BAMs: list[str], id: int, out:str, depth:int):
    outFile = os.path.join(out,i,".bam")    
    synProp = syn.getSyntheticLineageProportions(BAMs)
    synReads = syn.generateSyntheticReads(BAMs = BAMs, outFile = synBAMs, depth = depth)
    frejya = fn.runFrejya(BAMfile = synBAMs, variantOut = variantsOut, depthsOut = depthsOut, ref = ref, outFile = frejyaOut)
    os.remove(synReads)
    results = syn.compareSyntheticMix(frejya, synProp)
    return(results)

BAMs = fn.getBAMdb(seqPath = seqPath, metadata = metadata, BAMdb = BAMdb, BAMfiles = BAMfiles)
results = [createMix(BAMs, id, out, depth) for id in range(1,nMixes)]



