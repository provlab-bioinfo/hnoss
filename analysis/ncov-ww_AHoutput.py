import os, functions as fn
import pandas as pd
import searchTools as st
from pango_aliasor.aliasor import Aliasor
from datetime import datetime, date
os.chdir(os.path.dirname(__file__))

def generateAHoutput(freyjaFiles, IDs, output):
    # Merge Freyja and IDs
    ID_list = pd.read_excel(IDs)
    freyja = fn.formatFreyjaLineage(freyjaFiles).reset_index()
    freyja = ID_list.merge(freyja, on="file", how="left", indicator="Found")
    notFound = freyja.loc[freyja['Found'] == 'left_only']["SiteID"].values.tolist()
    print(f"Freyja output not found for ({len(notFound)}): {', '.join(notFound)}")
    freyja = freyja.drop(["Found"], axis=1, errors="ignore")
    # freyja = freyja.iloc[:, : 10]

    # Drop < 85% coverage and missing accessions
    lowCoverage = freyja[freyja.coverage < 85]["SiteID"].values.tolist()
    print(f"Strains dropped due to < 85% coverage ({len(lowCoverage)}): {', '.join(lowCoverage)}")
    freyja = freyja[freyja.coverage >= 85]
    freyja = freyja.dropna(subset=["Accession"])    

    # Add collection and reporting dates
    freyja["CollectionDate"] = freyja['SiteID'].apply(lambda x: datetime.strptime(x[3:9], '%y%m%d').date())
    freyja["AnalysisDate"] = freyja['Run'].apply(lambda x: datetime.strptime(x[0:6], '%y%m%d').date())
    freyja = freyja.drop(["Run","file","Notes","resid","coverage"], axis=1, errors="ignore")
    freyja.insert(2, "CollectionDate", freyja.pop("CollectionDate"))
    freyja.insert(3, "AnalysisDate", freyja.pop("AnalysisDate"))

    # Transform to long form
    # Not sure why doesn't work:
    # freyja = freyja.set_index(["Accession","SiteID"]).add_prefix("lineage-").reset_index()
    # pd.wide_to_long(freyja, stubnames = "lineage", i = ["Accession", "SiteID"], j = "strain", sep="-", suffix='\w+')
    idxCols = ["Accession","SiteID","CollectionDate","AnalysisDate"]
    freyja = freyja.set_index(idxCols)
    strainCols = freyja.columns
    freyja = freyja.reset_index()
    freyja = pd.melt(freyja, id_vars=idxCols, value_vars=strainCols)
    freyja = freyja.rename(columns={"variable": "Strain", "value": "Proportion"})
    freyja = freyja.dropna(subset=["Proportion"])

    # Add count and long_strain
    freyja['Count'] = freyja.groupby('SiteID')['SiteID'].transform('count')
    aliasor = Aliasor()
    freyja["long_strain"] = freyja['Strain'].apply(lambda x: aliasor.uncompress(x))

    # Export
    freyja = freyja.rename(columns={"Accession": "Sample"})
    freyja = freyja.sort_values(by=['CollectionDate',"SiteID"])
    freyja.to_csv(output, index = False)

ID_list = "/nfs/APL_Genomics/apps/development/ncov-wastewater/cumulative_ww_out.xlsx"
freyja = st.generateFlatFileDB("/nfs/APL_Genomics/apps/development/ncov-wastewater/freyja_agg_output/")
output = os.path.join("/nfs/APL_Genomics/apps/development/ncov-wastewater/", str(date.today().strftime("%Y%m%d")) + "_ww_out_test.csv")
generateAHoutput(freyjaFiles = freyja, IDs = ID_list, output = output)