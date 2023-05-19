from functions import *

os.chdir(os.path.dirname(__file__))
freyjaOut = ["./data/freyja1.tsv","./data/freyja2.tsv"]
fileCol = "file"
lineageCol = "lineages"
abundCol = "abundances"
lineages = pd.DataFrame()

def getLongForm(file):
    lineage = getFreyjaLineageProportions(file)
    lineage[fileCol] = os.path.basename(file)
    return(lineage)

lineages = pd.concat([getLongForm(file) for file in freyjaOut]).reset_index(drop=True)
lineages = lineages.groupby([fileCol, lineageCol])[abundCol].first().unstack()
cols = lineages.columns.tolist()
lineages["Dominant"] = lineages[cols].idxmax(axis=1)
lineages["Dominant%"] = lineages[cols].max(axis='columns')
lineages["SD"] = lineages[cols].std(axis=1)





print(lineages)

