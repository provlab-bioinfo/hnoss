import os, functions as fn
import pandas as pd
import searchTools as st
from pango_aliasor.aliasor import Aliasor
from datetime import datetime, date
from functools import reduce
os.chdir(os.path.dirname(__file__))

# WE accession
# Site ID
# Location
# Collection Date
# Sequencing Date
# Sample Key
# BAM file path
# Freyja variant file path
# Freyja demix file path
# Summarized Lineages
# Raw Lineages

def generateAHoutput(freyjaFiles, IDs, output):
    # Merge Freyja and IDs
    ID_list = pd.read_excel(IDs)
    freyjaRaw = fn.formatFreyjaLineage(freyjaFiles).add_suffix("_RawLineages")
    freyjaSumm = fn.formatFreyjaLineage(freyjaFiles, summarized=True).add_suffix("_SummarizedLineages")
    rawCols = freyjaRaw.columns.to_list()
    summCols = freyjaSumm.columns.to_list()
    freyja = freyjaSumm.join(freyjaRaw).reset_index()
    freyja = ID_list.merge(freyja, on="file", how="left", indicator="Found")
    notFound = freyja.loc[freyja['Found'] == 'left_only']["SiteID"].values.tolist()
    print(f"Freyja output not found for ({len(notFound)}): {', '.join(notFound)}")
    freyja = freyja.drop(["Found"], axis=1, errors="ignore")

    # Add collection and reporting dates
    freyja["CollectionDate"] = freyja['SiteID'].apply(lambda x: datetime.strptime(x[3:9], '%y%m%d').date())
    freyja["AnalysisDate"] = freyja['Run'].apply(lambda x: datetime.strptime(x[0:6], '%y%m%d').date())
    # freyja = freyja.drop(["Run","file","Notes","resid","coverage"], axis=1, errors="ignore")
    freyja.insert(3, "CollectionDate", freyja.pop("CollectionDate"))
    freyja.insert(4, "AnalysisDate", freyja.pop("AnalysisDate"))

    # Generate the Multi-index
    dataCols = list(set(freyja.columns) - set(rawCols) - set(summCols))
    freyja = freyja.rename(columns={c: c+'_SampleInfo' for c in freyja.columns if c in dataCols})
    reverseSplit = lambda str : '_'.join(str.split('_')[::-1])
    freyja.columns = [reverseSplit(col) for col in freyja.columns.to_list()]
    idx = freyja.columns.str.split('_', expand=True)
    freyja.columns = idx
    idx = pd.MultiIndex.from_product([idx.levels[0], idx.levels[1]])
    freyja.reindex(columns=idx, fill_value=-1)

    # Export
    freyja = freyja.sort_values(by=[('SampleInfo','Run'),("SampleInfo","CollectionDate")])
    freyja.to_csv(output, index = False)

ID_list = "/nfs/APL_Genomics/apps/development/ncov-wastewater/cumulative_ww_out.xlsx"
freyja = st.generateFlatFileDB("/nfs/APL_Genomics/apps/development/ncov-wastewater/freyja_agg_output/")
output = os.path.join("/nfs/APL_Genomics/apps/development/ncov-wastewater/", str(date.today().strftime("%Y%m%d")) + "_ww_out_wide.csv")
generateAHoutput(freyjaFiles = freyja, IDs = ID_list, output = output)