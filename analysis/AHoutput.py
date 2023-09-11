import os, functions as fn
import pandas as pd
from pango_aliasor.aliasor import Aliasor
os.chdir(os.path.dirname(__file__))

def generateAHoutput(freyjaFiles, IDs, output):
    # Merge Freyja and IDs
    ID_list = pd.read_excel(IDs)
    freyja = fn.formatFreyjaLineage(freyjaFiles)
    freyja = ID_list.merge(freyja, on="file")
    freyja = freyja.drop(["Run","file"], axis=1) #'resid', 'coverage',
    #freyja = freyja.iloc[:, : 10]

    # Transform to long form
    # Not sure why doesn't work:
    # freyja = freyja.set_index(["Accession","SiteID"]).add_prefix("lineage-").reset_index()
    # pd.wide_to_long(freyja, stubnames = "lineage", i = ["Accession", "SiteID"], j = "strain", sep="-", suffix='\w+')

    freyja = freyja.set_index(["Accession","SiteID"])
    cols = freyja.columns
    freyja = freyja.reset_index()
    freyja = pd.melt(freyja, id_vars=["Accession","SiteID"], value_vars=cols)
    freyja = freyja.rename(columns={"variable": "Strain", "value": "Proportion"})
    freyja = freyja.dropna()

    # Add count and long_strain
    freyja['Count'] = freyja.groupby('Accession')['Accession'].transform('count')
    aliasor = Aliasor()
    freyja["long_strain"] = freyja['Strain'].apply(lambda x: aliasor.uncompress(x))

    # Export
    freyja = freyja.rename(columns={"Accession": "Sample"})
    freyja = freyja.sort_values(by=['Sample'])
    freyja.to_csv(output, index = False)

ID_list = "/nfs/APL_Genomics/virus_covid19/ww_covid/cumulative output/cumulative_ww_out.xlsx"
freyja = "/nfs/APL_Genomics/virus_covid19/ww_covid/230907_N_I_004_WW/results/freyja/demix/230907_N_I_004_WW_freyja_agg.tsv"
output = "../results/230907_N_I_004_WW.csv"
generateAHoutput(freyjaFiles = freyja, IDs = ID_list, output = output)