import pandas as pd
import numpy as np
from pathlib import Path
import re
from functional import compose
from functools import reduce, partial
import itertools as it


# Speed check
# Check that numpy is compiled with OpenBLAS, Apple Accelerate or Intel MKL support
#np.__config__.show()

##################################################
#Helper High Order Functions
#map that return a list
def lmap(func, *iterable): return list(map(func, *iterable))

#compose list of functions
def mcompose(*func): return reduce(compose, func)

#################################################
#Part I: Functions for retrieving and dumping the data
#Composed functions / Functions dependant of data:

#Filter unused columns from Data Frame
def dfFilter(df, regexpKeepCol):
    return df.select(
            lambda colname:
            regexpKeepCol.search(colname), axis=1)

#1. Parse file path, extract (Experiment,Day,Cell line) from folder names
#3. Convert all values from string to integer
def xpparam(stringpath, regexpFolder):
    return compose(
        partial(lmap, int),
        regexpFolder.findall
    )(stringpath)

# 1. Read Excel/Csv
# 2. Filter out unused columns
# 3. Transform DataFrame to narrow format.
# Format is more suitable for automatic processing/Pivot tables and widely used for databases.
# Also format allows not to care about order of the data.
#    A  B            Data   Value
# 1  p  q    ==>  1    A     p
# 2  r  s         2    A     r
#                 1    B     q
#                 2    B     s
def xptonarrow(stringpath, regexpKeepCol):
    return mcompose(
        partial(pd.melt, id_vars="Name", var_name="Type"),
        partial(dfFilter, regexpKeepCol=reKeepCol),
        partial(pd.read_csv, delimiter=";")
        )(stringpath)

#Create tuple of ((Experiment, Day, Cell line), Data Frame(Std Deviation)
#in record/narrow format from list of Excel files
#Arg: list of paths, Regexp to extract metadata from path, Regexp to keep specific column names
def edcl(paths, reFpath, reFcol):
    return [(xpparam(str(path), reFpath), xptonarrow(str(path), reFcol)) for path in paths]

#Heavy lifting of the +/-/low/high/neg values
#Additional lifting: convert Replicat to Integer
def mapDico(*lconfig):
    if len(*lconfig) != 4:
        return None, None, None, None
    else:
        lconfig = tuple(*lconfig) #Can't use list directly as it's an unhashable type
        #CD24, CD44, Replicat, Comments/error
        return dicoRef[lconfig[0].lower()], dicoRef[lconfig[1].lower()], int(lconfig[2]), lconfig[3]

#Transform ( (a,b,c), DataFrame ) into a DataFrame with new columns A, B, C filled with values a,b,c
def mapMetaData(data):
    experiment = data[0][0]
    day = data[0][1]
    cell = data[0][2]
    DF = data[1]

    DF["Experiment"], DF["Day"], DF["Cell_Line"] = experiment, day, cell

    return DF

#Column "Name" holds CD24/CD44 +/-, Replicat ID
def extractMetaFromName(df, reParseName):
    DF = df #use a temp variable to avoid side effects

    DF["CD24"], DF["CD44"], DF["Replica"], DF["Comment_Error"] = \
    zip(*DF["Name"].map(compose(mapDico, reParseName.findall)))

    return DF

#######################################
# Functions for part II

#Add empty column
def add_none_col(df, i):
    df[i] = np.nan
    return df


#######################################
# Setup

#Regexp to parse info from folder path (probably working only on Linux, Mac and not on windows):
reFolder = re.compile(r'(?<=#0)\d+(?=/)|(?<=_D)\d+(?=/)|(?<=Result_SUM)\d+(?=\.csv)')
#Regexp to filter kept columns
reKeepCol = re.compile(r'^SD*|^Name')
#Regexp to parse info from column "Name"
reParseName = re.compile(r'((?<=CD24)[^C]+(?=CD44)|(?<=CD44)[^_]+(?=_)|(?<=_)\d|(?<=_\d).*(?=$))')

#Dictionary for consistent referencing (transform input to lower case for insensitive use)
dicoRef = {"+":"+", "high":"+", "med":"med", "low":"-", "-":"-", "neg":"-"}

#Order of the resulting columns
dicoSort = {'SDasinhAPC_data_neg_APC':1, 'SDasinhPE_data_neg_APC':2, 'SDasinhAPC_data_pos_APC':3, 'SDasinhPE_data_pos_APC':4}

#CSVseparator
csvsep = ";"

#######################################################
# Part I - Proc retrieve all the data and dump them

files_path = Path('dataset_only_transition').glob('**/*.csv')
data = edcl(files_path, reFolder, reKeepCol)
data_transf = pd.concat([mapMetaData(dataset) for dataset in data], ignore_index=True)
dtset = extractMetaFromName(data_transf, reParseName)

#Errors correction, potential side-effects due to in-place modification
#Error 1. med/high combo is actually high/high
dtset.loc[(dtset.CD24 == "med") & (dtset.CD44 == "+"), ["CD24", "Comment_Error"]] = \
    "+","CD24 changed from med to +"
#Error 2. tagged data named SUM149 with tag "_error SUM159" are actually for 159 cell line
dtset.loc[(dtset.Cell_Line == 149) & (dtset.Comment_Error == "_error SUM159"), ["Cell_Line", "Comment_Error"]] = \
    159, "_error SUM159. Corrected changed cell line from 149 to 159"

#Export the raw records
writer = pd.ExcelWriter("GraphPad_formatted_data.xlsx", engine='xlsxwriter')
dtset.to_excel(writer, "Raw data Std Dev in Asinh space")

#Retransform to linear space and export
dtlin = dtset
dtlin['value'] = dtlin['value'].map(np.sinh)
dtlin.to_excel(writer, "Geom-like sinh SD asinh")


#######################################################
# Part II - Pivot the data for copy/paste into GraphPad
#Clean data from control data rows, i.e. if CD24 or CD44 or Replicat are empty, it's not a regular experiment
clean_dt = dtset[dtset.CD24 != ""]

#Pivot the data and generate control data before reformatting
pivot = clean_dt.pivot_table(values = ["value"], \
    index = ["Cell_Line", "CD24", "CD44", "Day"], \
    columns = ["Type", "Experiment", "Replica"])
pivot.to_excel(writer, "Std Dev Asinh - Raw Pivot")
pivot.applymap(np.sinh).to_excel(writer, "Sinh SD Asinh - Raw Pivot")

#Data needs all 139 days even if no experiment -> rebuild the row indexes
fulldays_combo = it.product((149,159), ("+","-","med"), ("+","-","med"), range(1, 1+max(clean_dt['Day'])))

#Filter out (CD24+,CD44-) combination which doesn't exist
final_combo = [(a,b,c,d) for (a,b,c,d) in fulldays_combo if (b,c) not in [("+","-"),("-","med"),("+","med"),("med","-"),("med","+")]]

### Rebuild Indexes and export to csv
#Reindex rows
pivot_tosort = pivot.reindex(final_combo)

#Add empty columns 'Replica 7' after each Experiment 59 (GraphPad quirk)
cl = lmap(lambda x:x.tolist(), pivot_tosort.columns.levels)[:-2] + [[59],[7]]  #Build list of column level names + Experiment 59 + Replica 7
pivot_tosort = reduce(add_none_col, it.product(*cl), pivot_tosort) #Create all empty Experiment 59/Replica 7 columns

#Sort columns
pivot_final = pivot_tosort.reindex_axis( \
    sorted(pivot_tosort.columns, \
        key=lambda x: (dicoSort[x[1]], x[2], x[3])), axis=1)

#Export
pivot_final.to_excel(writer, "SD Asinh Formatted Final Pivot")
pivot_final.applymap(np.sinh).to_excel(writer, "Sinh SD Asinh - Formatted Final")

#Highlight the important table
writer.sheets["Sinh SD Asinh - Formatted Final"].set_tab_color('red')

writer.save()

# END
#######################################################
