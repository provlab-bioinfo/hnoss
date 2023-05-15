import time, os, re, pandas as pd, subprocess, random

# region: Functions

def makeFlatFileDB(dir: str, regex: str = None, fileExt: str = None, outFile: str = None, maxFilesRead:int = 100000000, excludeDirs: list[str] = [], verbose: bool = True):
    """Finds all files that fit a regex in a specified folder
    :param dir: Directory to search
    :param regex: The regex to search by
    :param outFile: The output file path
    :param maxFilesRead: The maximum number of files to read, defaults to 1000
    :param excludeDirs: List of folders to exclude
    :param verbose: Print progress messages?, defaults to True
    """    
    nFasta = nFiles = speed = 0
    out = []
    lastCheck = startTime = time.time()
    excludeDirs = "|".join(excludeDirs)

    if outFile is not None:
        out = open(outFile,'w')

    if (verbose): print("Searching for files...")

    for root, dirs, files in os.walk(dir, topdown=True):
        dirs.sort(reverse=True)
        if nFasta >= maxFilesRead: break
        for file in files:
            nFiles += 1
            match = file.endswith((fileExt)) if (regex is None) else re.match(regex, file)
            if match:
                if re.search(excludeDirs, root) != None: continue
                path = str(root) + "/" + str(file)
                if outFile is not None: 
                    out.write(path + "\n")
                else: 
                    out.append(path)
                nFasta += 1
            if (nFiles % 1000 == 0): 
                speed = str(round(1000/(time.time() - lastCheck)))
                lastCheck = time.time()
            if (verbose): print("   Parsed {} files and found {} FASTA files ({}/s)                ".format(nFiles,nFasta,speed), end="\r")

    if (verbose): print("   Parsed {} files and found {} FASTA files ({}s)                 ".format(nFiles,nFasta,str(round(time.time() - startTime,2))), end="\r")

    return (out if outFile is None else outFile)

    
def subsetFlatFileDB(inFile: str, outFile: str = None, includeTerms: list[str] = [], excludeTerms: list[str] = []):
    """Subsets a flat file database. 
    :param inFile: The original database path
    :param outFile: The path to save the subset database in
    :param includeTerms: Strings that paths must include. 
    :param excludeTerms: Strings that paths must not include. Suggest using: "Freebayes", "qc", "work"
    """
    if isinstance(excludeTerms, str): excludeTerms = [excludeTerms]
    if isinstance(includeTerms, str): includeTerms = [includeTerms]

    out = []
    if outFile is not None:
        out = open(outFile,'w')

    with open(inFile) as inDB:
        for line in inDB:
            inc = any(include in line for include in includeTerms) if len(includeTerms) else True
            exc = not any(exclude in line for exclude in excludeTerms) if len(excludeTerms) else True
            if (inc and exc): 
                if outFile is not None: 
                    out.write(line)
                else: 
                    out.append(str.strip(line))
    
    return (out if outFile is None else outFile)

def getBAMdb(seqPath, metadata, BAMdb, BAMfiles) -> pd.DataFrame:
    nanoBamRE = ".primertrimmed.rg.sorted.bam"
    illBamRE = ".mapped.primertrimmed.sorted.bam"

    if (not os.path.isfile(BAMdb)):
        makeFlatFileDB(seqPath, fileExt=((".bam")), outFile=BAMdb, excludeDirs=["work","tmp_bam","qc","test","troubleshooting"])

    if (not os.path.isfile(BAMfiles)):
        BAMs = subsetFlatFileDB(BAMdb, includeTerms=[nanoBamRE, illBamRE])
        BAMs = pd.DataFrame(BAMs, columns=["BAMPath"])
        BAMs['Key'] = BAMs['BAMPath'].transform(lambda x: os.path.basename(x).split('.')[0])
        BAMs = BAMs.drop_duplicates(subset=['Key'], keep='last')
        metadata = pd.read_csv(metadata, low_memory=False) 
        BAMs = BAMs.merge(metadata, on='Key', how='left')
        BAMs = BAMs[["Key", "BAMPath", "current_lineage"]].dropna()
        BAMs.to_csv(BAMfiles)
    else:
        BAMs = pd.read_csv(BAMfiles) 

    return BAMs

def getSyntheticList(BAMs:pd.DataFrame, n:int = 100) -> tuple[pd.DataFrame, pd.DataFrame]:
    BAMs = BAMs.sample(n = min(len(BAMs),n))
    return BAMs

def getLineageProportions(BAMs: pd.DataFrame):
    lineage = pd.DataFrame({'count' : BAMs.groupby("current_lineage").size()}).reset_index()
    lineage['prop'] = lineage['count'].transform(lambda x: round(100*x/len(lineage),2))
    return lineage

def generateSyntheticReads(BAMs: pd.DataFrame, workDir: str, outFile: str, depth:int = 10, seed:str = "seed"):
    BAMpaths = BAMs["BAMPath"].values.tolist()
    BAMpaths.sort()
    with open(outFile,mode="wb") as out:
        for BAM in BAMpaths:
            subsample = depth/len(BAMs)
            command = ["samtools", "view", "-s", str(subsample), BAM]
            subprocess.run(command, stdout=out)     
    return(outFile)

# endregion

os.chdir(os.path.dirname(__file__))
seqPath = "/nfs/APL_Genomics/virus_covid19/routine_seq/2023_01_Runs"
metadata = "/nfs/Genomics_DEV/projects/nextstrain/EBS-parser/results/allData.csv"
BAMdb = "./data/BAMdb.txt"
BAMfiles = "./data/BAMfiles.csv"
variantsOut = "./results/variants.out" 
depthsOut = "./results/depth.out"
ref = "./data/ncov.fasta"
workDir = "./work"
synBAMs = "./results/synBAMs.bam"
frejyaOut = "./results/freyja.tsv"

BAMs = getBAMdb(seqPath = seqPath, metadata = metadata, BAMdb = BAMdb, BAMfiles = BAMfiles)
synList = getSyntheticList(BAMs, n=100)
synProp = getLineageProportions(synList)
synReads = generateSyntheticReads(BAMs = synList, outFile = synBAMs, depth = 0.1)
