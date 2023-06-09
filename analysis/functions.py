import time, os, re, pandas as pd, subprocess, tempfile, numpy as np
from datetime import datetime
from searchTools import *
from pango_aliasor.aliasor import Aliasor

#region: pandas static vars
fileCol = "file"
lineageCol = "lineages"
summarizedCol = "summarized"
parentLineageCol = "parent_lineages"
abundCol = "abundances"
dominantCol = "Dominant"
dominantPercentCol = "Dominant %"
sdCol = "SD"
locationCol = "Location"
collectDateCol = "Collection Date"
residualCol = "resid"
coverageCol = "coverage"
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

def aggregateFreyja(freyja, freyjaOut):    
    #freyja aggregate [directory-of-output-files] --output [aggregated-filename.tsv] --ext output
    subprocess.run(["freyja","aggregate",freyja,"--output", freyjaOut,"--ext","freyja.tsv"])
    return freyjaOut

def plotFreyja(freyja, output):
    #freyja plot [aggregated-filename-tsv] --output [plot-filename(.pdf,.png,etc.)]
    subprocess.run(["freyja","plot",freyja,"--output", "summ-" + output])
    subprocess.run(["freyja","plot",freyja,"--lineages","--output", "lineage-" + output])        

def startFreyjaDashboard(freyja, metadata, output):
    #freyja dash [aggregated-filename-tsv] [sample-metadata.csv] [dashboard-title.txt] [introContent.txt] --output [outputname.html]
    subprocess.run(["freyja","dash",freyja,metadata,"--output",output])

def formatFreyjaOutput(file:list[str]) -> pd.DataFrame:
    """Gets the lineage proportions from Frejya output files(s)
    :param file: The path to the Freyja output file(s)
    :return: A DataFrame with columns for lineages and abundances
    """    
    if (isinstance(file, list)): # For case of inputting individual files
        freyja = [pd.read_csv(f, sep="\t", index_col = 0) for f in file]
        freyja = [f.transpose() for f in freyja]
        freyja = pd.concat([freyja])
    else:
        freyja = pd.read_csv(file, sep="\t", index_col = 0)
        if (len(freyja.columns) == 1): # For a single file inputted
            freyja = freyja.transpose()

    freyja = freyja.rename_axis('file').reset_index()
    freyja[lineageCol] = freyja[lineageCol].str.split(" ")
    freyja[abundCol] = freyja[abundCol].str.split(" ")
    freyja = freyja.explode([lineageCol,abundCol]).reset_index(drop=True)
    freyja[abundCol] = freyja[abundCol].transform(lambda x: sigfig(100*float(x)))
    freyja = freyja.groupby([fileCol, residualCol, coverageCol, lineageCol])[abundCol].first().unstack()
    return(freyja)

def collapseFreyjaLineage(freyja: pd.DataFrame, strains: list[str]):
    """ Collapses Freyja lineages. Cannot parse down past "A","B" or any recombinant strains
    :param freyja: A Freyja output dataset formatted with formatFreyjaOutput()
    :param strains: The strains to parse down until
    """    
    aliasor = Aliasor()

    while(True):
        cols = list(freyja.columns)
        cols = [col for col in cols if col not in [fileCol,residualCol,coverageCol]]    
        badcols = set([strain for strain in cols if strain not in strains])
        beforeCols = set(freyja.columns)

        for strain in badcols: 
            parent = aliasor.parent(strain)
            if (parent == ""): continue
            freyja = freyja.rename(columns={strain: parent})
            
        freyja = freyja.groupby(freyja.columns, axis=1).sum()
        freyja = freyja.replace({0:np.nan})

        if (beforeCols == set(freyja.columns)): break

    return freyja

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
