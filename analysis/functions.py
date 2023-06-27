import time, os, re, pandas as pd, subprocess, tempfile, numpy as np
import configSettings as cfg
import searchTools as st
from datetime import datetime

from pango_aliasor.aliasor import Aliasor

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
        st.generateFlatFileDB(seqPath, fileExt=((".bam")), outFile=BAMdb)#, excludeDirs=["work","tmp_bam","qc","test","troubleshooting"])

    if (not os.path.isfile(BAMdata)):
        BAMs = st.searchFlatFileDB(BAMdb, includeTerms=[nanoBamRE, illBamRE])#,  excludeTerms=["work","ncovIllumina","bai","test","results_old","results_unmerged"])
        BAMs = pd.DataFrame(BAMs, columns=[cfg.BAMPathCol])       
        BAMs['Key'] = BAMs[cfg.BAMPathCol].transform(lambda x: os.path.basename(x).split('.')[0])
        BAMs = BAMs.drop_duplicates(subset=['Key'], keep='last')
        if (metadata is not None):
            metadata = pd.read_csv(metadata, low_memory=False) 
            BAMs = BAMs.merge(metadata, on='Key', how='left')
            BAMs = BAMs[["Key", cfg.BAMPathCol, "current_lineage","collection_date"]].dropna()
        else:
            BAMs = BAMs[["Key", cfg.BAMPathCol]].dropna()
        BAMs.to_csv(BAMdata, index=False)
    else:
        BAMs = pd.read_csv(BAMdata, index_col=False) 

    return BAMs
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

def formatFreyjaLineage(file:list[str]) -> pd.DataFrame:
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
    freyja[cfg.lineageCol] = freyja[cfg.lineageCol].str.split(" ")
    freyja[cfg.abundCol] = freyja[cfg.abundCol].str.split(" ")
    freyja = freyja.explode([cfg.lineageCol,cfg.abundCol]).reset_index(drop=True)
    freyja[cfg.abundCol] = freyja[cfg.abundCol].transform(lambda x: float(sigfig(100*float(x))))
    freyja = freyja.groupby([cfg.fileCol, cfg.residualCol, cfg.coverageCol, cfg.lineageCol])[cfg.abundCol].first().unstack()
    return(freyja)

def collapseFreyjaLineage(freyja: pd.DataFrame, strains: list[str]):
    """ Collapses Freyja lineages. Cannot parse down past "A","B" or any recombinant strains
    :param freyja: A Freyja output dataset formatted with formatFreyjaOutput()
    :param strains: The strains to parse down until
    """    
    aliasor = Aliasor()

    while(True):
        cols = list(freyja.columns)
        cols = [col for col in cols if col not in [cfg.fileCol, cfg.residualCol, cfg.coverageCol]]    
        badcols = set([strain for strain in cols if strain not in strains])
        beforeCols = set(freyja.columns)

        for strain in badcols: 
            freyja = collapseStrain(freyja, strain, aliasor)
            
        freyja = freyja.groupby(freyja.columns, axis=1).sum()
        freyja = freyja.replace({0:np.nan})

        if (beforeCols == set(freyja.columns)): break

    return freyja

def collapseStrain(freyja: pd.DataFrame, strain: str, aliasor:Aliasor = None):
    """Collapse a single strain in a formatted Freyja dataset
    :param freyja: A Freyja output dataset formatted with formatFreyjaOutput()
    :param strain: The name of the strain to collapse. If it's a terminal parent (e.g., 'A','B', any 'X' strain), it cannot be collapsed further.
    :param aliasor: A pango_aliasor object. Will be generated if not specified (adds significant runtime if done repeatedly), defaults to None
    :return: A Freyja dataset with the collapsed strain (if possible)
    """        
    if (aliasor == None): aliasor = Aliasor()
    parent = aliasor.parent(strain)
    if (parent != ""): freyja = freyja.rename(columns={strain: parent})
    return freyja

def collapseSamples(freyja:pd.DataFrame, date = True, location = True):
    """Collapses samples by dates and/or locations. Residuals and coverage are lost.
    :param freyja: A Freyja output dataset formatted with formatFreyjaOutput()
    :param date: Collapse by date?, defaults to True
    :param location: Collapse by location?, defaults to True
    :return: A DataFrame containing mean of each strain for dates and/or locations
    """    
    freyja.drop(columns=[cfg.coverageCol, cfg.residualCol])

    if (date and not location):
        freyja = freyja.groupby([cfg.collectDateCol]).mean()
    if (location and not date):
        freyja = freyja.groupby([cfg.locationCol]).mean()
    if (date and location):
        freyja = freyja.groupby([cfg.collectDateCol, cfg.locationCol]).mean()
    if (not date and not location):
        freyja = freyja.mean()

    return freyja

def locateFreyjaSamples(freyja: pd.DataFrame) -> pd.DataFrame: 
    """Adds location and collection date to Freyja DataFrame
    :param freyja: A DataFrame from collateFreyjaSamples()
    :return: A DataFrame with location and colletion data
    """    
    freyja = freyja.reset_index()
    freyja[cfg.locationCol] = freyja[cfg.fileCol].transform(lambda x: x[0:2])
    freyja[cfg.collectDateCol] = freyja[cfg.fileCol].transform(lambda x: datetime.strptime(x[3:9], '%y%m%d').strftime("%Y-%m-%d"))
    freyja = freyja.set_index(cfg.fileCol)
    return (freyja)
#endregion

def plotFreyja(freyja: pd.DataFrame) -> pd.DataFrame:
    freyja = collapseSamples(freyja, date = True, location = True)

def generateAuspiceFreqs(freyja:pd.DataFrame, outFile: str):
    """Generates tip_frequences.json for Auspice. Each entry will have frequency proportion only for the day the reading was taken
    :param freyja: a formatted Freyja DataFrame
    :param outFile: The output JSON
    """    
    freyja[cfg.weekCol] = freyja[cfg.collectDateCol].transform(lambda x: fn.dateToFractionalWeek(x))

    # Generate frequencies JSON
    freqs = sorted(list(set(freyja[cfg.weekCol].values)))
    freqs = pd.DataFrame(index = freyja.index.values.tolist(),columns = freqs)
    for index, row in freyja.iterrows(): freqs.at[index, row[cfg.weekCol]] = 1
    freqs.update(freqs.div(freqs.sum(axis=0),axis=1).fillna(0))
    freqs[cfg.freqCol] = freqs.values.tolist()
    pivots = list(freqs.columns)
    pivots.remove(cfg.freqCol)
    freqs = freqs[[cfg.freqCol]]
    json = loads(freqs.to_json(orient="index"))
    json.update(  {"generated_by": {
        "program": "custom",
        "version": "0.0.1"
    }})
    json.update({"pivots": pivots})
    json = dumps(json, indent=4)

    with open(outFile, "w") as out:
        out.write(json)

    # # Generate fake tree JSON
    # freyja["div"] = 0.01
    # freyja["num_date"] = freyja[cfg.weekCol].transform(lambda x: {"value": x,"confidence":[x,x]})
    # tree = {}

    # for index,row in freyja.iterrows():
    #     row = row.to_dict()
    #     print(row)
    #     exit()

#endregion

#region: accFuncs
def sigfig(val, n:int = 3):
    """Forces value to specific number of decimal points
    :param val: The value to format
    :param n: The number of decimal places
    :return: The truncated float
    """    
    # if (n == 0): return round(val)
    return '{0:.{1}f}'.format(float(val),n)
#endregion
