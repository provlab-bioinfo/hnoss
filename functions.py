import time, os, re, pandas as pd, subprocess, tempfile
from datetime import datetime
from searchTools import *

#region: pandas static vars
fileCol = "file"
lineageCol = "lineages"
parentLineageCol = "parent_lineages"
abundCol = "abundances"
dominantCol = "Dominant"
dominantPercentCol = "Dominant %"
sdCol = "SD"
locationCol = "Location"
collectDateCol = "Collection Date"
residualCol = "Residual"
coverageCol = "Coverage"
BAMPathCol = "BAMPath"
#endregion

#region: database
def getBAMdb(seqPath: str, BAMdb, BAMdata, metadata: str = None) -> pd.DataFrame:
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
        generateFlatFileDB(seqPath, fileExt=((".bam")), outFile=BAMdb)#, excludeDirs=["work","tmp_bam","qc","test","troubleshooting"])

    if (not os.path.isfile(BAMdata)):
        BAMs = searchFlatFileDB(BAMdb, includeTerms=[nanoBamRE, illBamRE])#,  excludeTerms=["work","ncovIllumina","bai","test","results_old","results_unmerged"])
        BAMs = pd.DataFrame(BAMs, columns=[BAMPathCol])       
        BAMs['Key'] = BAMs[BAMPathCol].transform(lambda x: os.path.basename(x).split('.')[0])
        BAMs = BAMs.drop_duplicates(subset=['Key'], keep='last')
        if (metadata is not None):
            metadata = pd.read_csv(metadata, low_memory=False) 
            BAMs = BAMs.merge(metadata, on='Key', how='left')
            BAMs = BAMs[["Key", BAMPathCol, "current_lineage","collection_date"]].dropna()
        else:
            BAMs = BAMs[["Key", BAMPathCol]].dropna()
        BAMs.to_csv(BAMdata, index=False)
    else:
        BAMs = pd.read_csv(BAMdata, index_col=False) 

    return BAMs
#endregion

#region: synthetic samples
def getSyntheticList(BAMs:pd.DataFrame, n:int = 100) ->  pd.DataFrame:
    """Creates a random sample of COVID BAM entries
    :param BAMs: The dataframe containing the BAM information. Contains columns 'Key', 'BAMpath', 'current_lineage'
    :param n: the number of BAM files to randomly select, defaults to 100
    :return: The subsetted dataframe
    """    
    BAMs = BAMs.sample(n = min(len(BAMs),n))
    return BAMs

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

def getSyntheticLineageProportions(BAMs: pd.DataFrame, parentLineage:bool = True) -> pd.DataFrame:
    """Gets a summary of lineage proportions for COVID BAM samples
    :param BAMs: The dataframe containing the BAM information. Contains columns 'Key', 'BAMpath', 'current_lineage'
    :return: A dataframe containing columns for 'lineages' and 'abundances' of COVID strains
    """    
    lineage = pd.DataFrame({'count' : BAMs.groupby("current_lineage").size()}).reset_index()
    lineage = lineage.rename(columns={"current_lineage": lineageCol})
    samples = lineage['count'].sum()
    lineage[abundCol] = lineage['count'].transform(lambda x: sigfig(100*(x/samples)))
    lineage = lineage.drop(columns=['count'])
    if (parentLineage): lineage = getParentLineage(lineage)
    return lineage

def compareFreyja(freyjaOut, lineage):
    """Compares a Freyja analysis to a synthetic input
    :param freyjaOut: The path to the Frejya output TSV
    :param lineage: The lineage of a synthetic sample, from getLineageProportions()
    :return: A dataframe containing the abundances (%) of each variant and parent variant 
    """    
    freyja = getFreyjaLineageProportions(freyjaOut)
    freyja = freyja.merge(lineage, how="outer", on=lineageCol, suffixes=["_freyja","_synthetic"])
    return(freyja)
#endregion

#region: Freyja
def runFrejya(BAMfile: list[str], outDir: str, ref: str, refname: str = None) -> str:
    """Runs a Freyja analysis on mixed COVID samples. See https://github.com/andersen-lab/Freyja
    :param BAMfile: The BAM files to de-mix    
    :param outDir: The output directory
    :param ref: The genome FASTA reference
    """    
    if isinstance(BAMfile, list):
        freyjaOut = [runFrejya(file, outDir, ref, refname) for file in BAMfile]
    else:
        getOutFile = lambda ext: os.path.join(outDir, os.path.splitext(os.path.basename(BAMfile))[0] + ext)
        variantsOut = getOutFile(".variants.tsv")
        depthsOut = getOutFile(".depths.tsv")
        freyjaOut = getOutFile(".freyja.tsv")

        # freyja variants [bamfile] --variants [variant outfile name] --depths [depths outfile name] --ref [reference.fa]
        cmd = ["freyja", "variants", BAMfile, "--variants", variantsOut, "--depths", depthsOut, "--ref", ref]
        #if (refname is not None): cmd = cmd + ["--refname", refname]
        subprocess.run(cmd)

        # Not sure why necessary. Cannot use --refname in above command, or index file is not found, and cannot have zero count reads, or below command has a log error
        depths = pd.read_csv(depthsOut, names = ["ref","pos","nt","count"], sep="\t")
        depths = depths[depths.ref.isin([refname])]
        depths.to_csv(depthsOut, header=None, index=None, sep="\t")

        # freyja demix [variants-file] [depth-file] --output [output-file]
        subprocess.run(["freyja","demix",variantsOut,depthsOut,"--output",freyjaOut])
    return freyjaOut    

def getParentLineage(abundances: pd.DataFrame) -> pd.DataFrame:
    """Get parent strain lineages from Freyja output
    :param abundances: A DataFrame containing full lineages and abundances
    :return: A DataFrame containing only parent lineages and abundances
    """    
    abundances[parentLineageCol] = abundances[lineageCol].transform(lambda x: "VF_" + x.removeprefix("Var_").split('.')[0])
    abundances = abundances.drop(columns=[lineageCol]).groupby(parentLineageCol).sum().reset_index()
    abundances = abundances.sort_values(by=[abundCol]).rename(columns={parentLineageCol: lineageCol})
    return(abundances)

def getFreyjaLineageProportions(file:str, parentLineage:bool = True) -> pd.DataFrame:
    """Gets the lineage proportions from a Freyja output file
    :param file: The path to the Freyja output file
    :param parentLineage: Should parent lineages also be included?, defaults to True
    :return: A DataFrame with columns for lineages and abundances
    """    
    freyja = pd.read_csv(file, sep="\t", index_col=0)
    freyja = freyja.loc[[lineageCol,abundCol]].transpose()
    freyja[lineageCol] = freyja[lineageCol].str.split(" ")
    freyja[abundCol] = freyja[abundCol].str.split(" ")
    freyja = freyja.explode([lineageCol,abundCol]).reset_index(drop=True)
    freyja[abundCol] = freyja[abundCol].transform(lambda x: sigfig(100*float(x)))
    freyja[lineageCol] = freyja[lineageCol].transform(lambda x: "Var_" + x)
    if (parentLineage): freyja = pd.concat([freyja, getParentLineage(freyja)])
    return(freyja)

def getFreyjaConfidence(freyjaPath:str) -> pd.DataFrame:
    """Gets residual and coverage from Freyja output files
    :param file: The path to the Freyja output file
    :return: A DataFrame with columns for residuals and coverage
    """    
    freyja = [pd.read_csv(file, sep="\t", index_col=0) for file in freyjaPath]
    freyja = pd.concat([f.loc[["resid","coverage"]].transpose() for f in freyja])
    freyja = freyja.rename(columns={"resid": residualCol, "coverage": coverageCol})
    freyja = freyja.applymap(sigfig)
    freyja[fileCol] = [os.path.basename(path) for path in freyjaPath]
    freyja = freyja.set_index(fileCol)
    return (freyja)

def collateFreyjaSamples(freyjaPath: list[str]) -> pd.DataFrame:
    """Collates wastewater samples into long form table
    :param freyjaOut: Paths to the Freyja output files
    """    
    def getLongForm(file):
        lineage = getFreyjaLineageProportions(file)
        lineage[fileCol] = os.path.basename(file)
        return(lineage)

    freyja = pd.concat([getLongForm(file) for file in freyjaPath]).reset_index(drop=True)
    freyja = freyja.groupby([fileCol, lineageCol])[abundCol].first().unstack()
    return (freyja)

def getFreyjaVariantStats(freyja: pd.DataFrame) -> pd.DataFrame:
    """Adds lineage statistics to a long form Freyja table
    :param freyja: A DataFrame from collateFreyjaSamples()
    :return: A DataFrame with variant statistics
    """    
    cols = freyja.columns.tolist()
    cols = [i for i in cols if "Var_" in i]
    freyja[dominantCol] = freyja[cols].idxmax(axis=1)
    freyja[dominantPercentCol] = freyja[cols].max(axis='columns')
    freyja[sdCol] = freyja[cols].std(axis=1)
    freyja[sdCol] = freyja[sdCol].apply(sigfig)
    return (freyja)

def locateFreyjaSamples(freyja: pd.DataFrame) -> pd.DataFrame: 
    """Adds location and collection date to Freyja DataFrame
    :param freyja: A DataFrame from collateFreyjaSamples()
    :return: A DataFrame with location and colletion data
    """    
    freyja = freyja.reset_index()
    freyja[locationCol] = freyja[fileCol].transform(lambda x: x[0:2])
    freyja[collectDateCol] = freyja[fileCol].transform(lambda x: datetime.strptime(x[3:9], '%y%m%d').strftime("%Y-%m-%d"))
    freyja = freyja.set_index(fileCol)
    return (freyja)

def formatFreyjaOutput(freyjaPath):
    freyja = collateFreyjaSamples(freyjaPath)
    freyjaConf = getFreyjaConfidence(freyjaPath)
    freyja = freyja.merge(freyjaConf, on=fileCol)
    freyja = getFreyjaVariantStats(freyja)
    #freyja = locateFreyjaSamples(freyja)
    return(freyja)
#endregion

#region: accFuncs
def sigfig(val, n:int = 3):
    """Forces value to specific number of decimal points
    :param val: The value to format
    :param n: The number of decimal places
    :return: The truncated float
    """    
    return float('{0:.{1}f}'.format(float(val),n))
#endregion