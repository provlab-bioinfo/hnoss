import time, os, re, pandas as pd, subprocess, random, tempfile

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

def getBAMdb(seqPath: str, metadata, BAMdb, BAMfiles) -> pd.DataFrame:
    """Generates list of BAM files for a specific folder.
    Uses regular expressions for the expected BAM output files from Nanopore and Illumina sequencing
    :param seqPath: The path to the folder containing BAM files
    :param metadata: The path to the CSV containing patient metadata
    :param BAMdb: The path to the output database. Will generate it if the file is not present.
    :param BAMfiles: The path to output the collated BAM file and patient metadata. Will generate it if the file is not present.
    :return: A dataframe containing 'Key', 'BAMpath', 'current_lineage'
    """    
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

def getSyntheticList(BAMs:pd.DataFrame, n:int = 100) ->  pd.DataFrame:
    """Creates a random sample of COVID BAM entries
    :param BAMs: The dataframe containing the BAM information. Contains columns 'Key', 'BAMpath', 'current_lineage'
    :param n: the number of BAM files to randomly select, defaults to 100
    :return: The subsetted dataframe
    """    
    BAMs = BAMs.sample(n = min(len(BAMs),n))
    return BAMs

def getLineageProportions(BAMs: pd.DataFrame) -> pd.DataFrame:
    """Gets a summary of lineage proportions for COVID BAM samples
    :param BAMs: The dataframe containing the BAM information. Contains columns 'Key', 'BAMpath', 'current_lineage'
    :return: A dataframe containing columns for 'lineages' and 'abundances' of COVID strains
    """    
    lineage = pd.DataFrame({'count' : BAMs.groupby("current_lineage").size()}).reset_index()
    lineage = lineage.rename(columns={"current_lineage": "lineages"})
    samples = lineage['count'].sum()
    lineage['abundances'] = lineage['count'].transform(lambda x: sigfig(100*(x/samples)))
    return lineage.drop(columns=['count'])

def generateSyntheticReads(BAMs: pd.DataFrame, outFile: str, depth:int = 10, seed:str = "seed") -> str:
    """Generates a synthetic BAM file from a list of COVID samples
    :param BAMs: The dataframe containing the BAM information. Contains columns 'Key', 'BAMpath', 'current_lineage'
    :param outFile: The path to the output BAM file
    :param depth: The depth of sampling for the input files. For each sample [depth]/[# of samples] = % of total reads , defaults to 10
    :return: outFile; the path to the output BAM file
    """    
    BAMpaths = BAMs["BAMPath"].values.tolist()
    BAMpaths.sort()

    with tempfile.NamedTemporaryFile() as tmp, open(outFile,mode="wb") as out:
        for idx, BAM in enumerate(BAMpaths):
            subsample = depth/len(BAMs)
            if (idx == 0): command = ["samtools", "view", "-h","-s", str(subsample), BAM]
            else: command = ["samtools", "view", "-s", str(subsample), BAM]
            subprocess.run(command, stdout=tmp)     
        command = ["samtools", "sort", tmp.name]
        subprocess.run(command, stdout=out)   
        tmp.close()    
    
    return(outFile)

def runFrejya(BAMfile: str, variantOut: str, depthsOut: str, ref: str, outFile:str):
    """Runs a Freyja analysis on mixed COVID samples. See https://github.com/andersen-lab/Freyja
    :param BAMfile: The BAM file to de-mix    
    :param variantOut: The path to the output variant file
    :param depthsOut: The path to the output depths file
    :param outFile: The path to the output Freyja file
    :return: outFile: The path to the output Freyja file 
    """    
    # freyja variants [bamfile] --variants [variant outfile name] --depths [depths outfile name] --ref [reference.fa]
    subprocess.run(["freyja", "variants", BAMfile, "--variants", variantOut, "--depths", depthsOut, "--ref", ref])

    # freyja demix [variants-file] [depth-file] --output [output-file]
    subprocess.run(["freyja","demix",variantOut,depthsOut,"--output",outFile])

    return outFile

def compareFreyja(freyjaOut, lineage):
    """Compares a Freyja analysis to a synthetic input
    :param freyjaOut: The path to the Frejya output TSV
    :param lineage: The lineage of a synthetic sample, from getLineageProportions()
    :return: A dataframe containing the abundances (%) of each variant and parent variant 
    """    
    freyja = pd.read_csv(freyjaOut, sep="\t", index_col=0)
    freyja = freyja.loc[["lineages","abundances"]].transpose()
    freyja["lineages"] = freyja["lineages"].str.split(" ")
    freyja["abundances"] = freyja["abundances"].str.split(" ")
    freyja = freyja.explode(["lineages","abundances"]).reset_index(drop=True)
    freyja['abundances'] = freyja['abundances'].transform(lambda x: sigfig(100*float(x)))
    freyja = freyja.merge(lineage, how="outer", on="lineages", suffixes=["_freyja","_synthetic"])
    freyja['parent_lineages'] = freyja['lineages'].transform(lambda x: x.split('.')[0])
    freyja2 = freyja.groupby('parent_lineages').sum().reset_index()
    freyja2.rename(columns={'lineages':'lineagesBAD', 'parent_lineages':'lineages'}, inplace=True)
    freyja2.rename(columns={'lineagesBAD':'parent_lineages'}, inplace=True)
    freyja2.sort_values(by=['abundances_freyja'], inplace=True)
    freyja = pd.concat([freyja, freyja2]).drop(columns=['parent_lineages'])
    return(freyja)

def sigfig(val, n:int = 3):
    """Forces value to specific number of decimal points
    :param val: The value to format
    :param n: The number of decimal places
    :return: The truncated float
    """    
    return float('{0:.{1}f}'.format(float(val),n))

# endregion

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
synProp = getLineageProportions(synList)
synReads = generateSyntheticReads(BAMs = synList, outFile = synBAMs, depth = 10)
frejya = runFrejya(BAMfile = synBAMs, variantOut = variantsOut, depthsOut = depthsOut, ref = ref, outFile = frejyaOut)

#synProp = pd.read_csv(synPropOut, sep="\t")
results = compareFreyja(frejyaOut, synProp)
print(results)
results.to_csv(compareOut, sep="\t", index = False)