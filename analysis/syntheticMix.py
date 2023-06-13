import functions as fn
import configSettings as cfg
import pandas as pd, tempfile, subprocess

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

def getSyntheticLineageProportions(BAMs: pd.DataFrame) -> pd.DataFrame:
    """Gets a summary of lineage proportions for COVID BAM samples
    :param BAMs: The dataframe containing the BAM information. Contains columns 'Key', 'BAMpath', 'current_lineage'
    :return: A dataframe containing columns for 'lineages' and 'abundances' of COVID strains
    """    
    lineage = pd.DataFrame({'count' : BAMs.groupby("current_lineage").size()}).reset_index()
    lineage = lineage.rename(columns={"current_lineage": cfg.lineageCol})
    samples = lineage['count'].sum()
    lineage[cfg.abundCol] = lineage['count'].transform(lambda x: float(fn.sigfig(100*(x/samples))))
    lineage = lineage.drop(columns=['count'])
    return lineage

def compareSyntheticMix(freyja, syntheticMix):
    """Compares a Freyja analysis to a synthetic input
    :param freyjaOut: The path to the Frejya output TSV
    :param lineage: The lineage of a synthetic sample, from getLineageProportions()
    :return: A dataframe containing the abundances (%) of each variant and parent variant 
    """    
    freyja = fn.formatFreyjaLineage(freyja)
    freyja = fn.collapseFreyjaLineage(freyja, syntheticMix["current_lineage"].value.tolist())
    freyja = freyja.merge(syntheticMix, how="outer", on=cfg.lineageCol, suffixes=["_freyja","_synthetic"])
    return(freyja)
