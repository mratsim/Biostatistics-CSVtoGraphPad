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


#Lifting: convert to integer
def mapDico(*lconfig):
    if len(*lconfig) != 2:
        return None, None
    else:
        lconfig = tuple(*lconfig) #Can't use list directly as it's an unhashable type
        #IniCells, Replica
        return int(lconfig[0]), int(lconfig[1])

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

    DF["IniCells"], DF["Replica"] = \
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
reParseName = re.compile(r'(?<=[T2]_)\d+(?=_)|(?<=0_)\d(?=$)')

#Order of the resulting columns
dicoSort = {'SDasinhAPC_data_neg_APC':1, 'SDasinhPE_data_neg_APC':2, 'SDasinhAPC_data_pos_APC':3, 'SDasinhPE_data_pos_APC':4}

#CSVseparator
csvsep = ";"

#######################################################
# Part I - Proc retrieve all the data and dump them

files_path = Path('dataset_only_tumorsphere').glob('**/*.csv')
data = edcl(files_path, reFolder, reKeepCol)
data_transf = pd.concat([mapMetaData(dataset) for dataset in data], ignore_index=True)
dtset = extractMetaFromName(data_transf, reParseName)

#Export the raw records
writer = pd.ExcelWriter("GraphPad_fmt_data_tumorsphere_v2.xlsx", engine='xlsxwriter')
dtset.to_excel(writer, "Raw data Std Dev in Asinh space")

#Retransform to linear space and export
dtlin = dtset
dtlin['value'] = dtlin['value'].map(np.sinh)
dtlin.to_excel(writer, "Geom-like sinh SD asinh")


#######################################################
# Part II - Pivot the data for copy/paste into GraphPad
#Clean data from control data rows, i.e. if CD24 or CD44 or Replicat are empty, it's not a regular experiment
clean_dt = dtset[dtset.IniCells != ""]

#Pivot the data and generate control data before reformatting
pivot = clean_dt.pivot_table(values = ["value"], \
    index = ["Cell_Line", "IniCells", "Day"], \
    columns = ["Type", "Experiment", "Replica"])
pivot.to_excel(writer, "Std Dev Asinh - Raw Pivot")
pivot.applymap(np.sinh).to_excel(writer, "Sinh SD Asinh - Raw Pivot")

#Add empty columns 'Replica 5' after each Experiment 91 (GraphPad quirk)
cl = lmap(lambda x:x.tolist(), pivot.columns.levels)[:-2] + [[91],[5]]  #Build list of column level names + Experiment 59 + Replica 7
pivot_tosort = reduce(add_none_col, it.product(*cl), pivot) #Create all empty Experiment 59/Replica 7 columns

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
