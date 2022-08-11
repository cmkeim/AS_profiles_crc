# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 22:20:51 2022

@author: cmk_8
"""

# -*- coding: utf-8 -*-
"""
Created on Sat May 21 12:04:11 2022

@author: cmk_8
"""
# https://crema.unibas.ch/CRUMARA/jobs/data_uokaj6nh/report/status.html
import pandas as pd

# Get .psi data from TCGA SPliceSeq web application (psi_reader.py)
read_psi = pd.read_csv('C:/Users/cmk_8/OneDrive/data/PSI_download_READ.txt', sep="\t")
read_psi["cancer_type"] = "READ"
coad_psi = pd.read_csv("C:/Users/cmk_8/OneDrive/data/PSI_download_COAD.txt", sep="\t")
coad_psi["cancer_type"] = "COAD"

crc_psi = coad_psi.append(read_psi)
# crc_psi.to_csv("C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/PSI_COAD_and_READ.tsv", sep="\t", index=False)


########################################################################################
# Import full set of STRs -- 17.6.2022
########################################################################################
str_full = pd.read_csv('C:/Users/cmk_8/OneDrive/data/tral_and_perf_panel.tsv', header=None, sep="\t", names=[
    "chromosome", "start", "end", "repeat_length", "msa", "strand", "db_match"])
str_full.head()

# Change/Remove order of columns to receive a .bed file --> use for liftOver to hg19 coordinates
# new_columns_str = ["chromosome", "start", "end", "msa","repeat_length", "strand"]
# str_full_bed = str_full[new_columns_str]
# str_full_bed.to_csv("C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/str_full_bed.bed", 
#                                sep="\t", header=False, index=False)

# -------- Import the full "str_perf_panel" liftover to hg19-coordinates    
# --------  NEW 27.6.22 - STR panel from hg38 to hg19 lifted coordinates
# liftOver -bedPlus=3 /mnt/c/Users/cmk_8/OneDrive/data/str_full_bed2.bed  /mnt/c/Users/cmk_8/OneDrive/data/hg38ToHg19.over.chain.gz
str_full_lift = pd.read_csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/str_full_Hg38ToHg19.bed', header=None, sep="\t", names=[
    "chromosome", "start", "end", "msa", "repeat_length", "strand", "db_match"])

    
# ########################################################################################
# # NEW FILTERING - Filter the liftOver with the original .psi data
# ########################################################################################

# This will filter the dataframe according the matching gene symbol
df_hg19_filtered_crc = df_hg19[df_hg19.symbol.isin(crc_psi.symbol)]

# This will filter the dataframe according the matching gene symbol
df_hg38_lift_filtered_crc = df_annotation_lift[df_annotation_lift.symbol.isin(crc_psi.symbol)]

# CONTROL: should match the crc_unique
# df_hg19_filtered_crc = df_hg19_filtered_crc.Symbol.unique()
# df_hg38_lift_filtered_crc = df_hg38_lift_filtered_crc.symbol.unique()
# crc_unique = sorted(crc_psi.symbol.unique())

########################################################################################
# Create .bed file for bedtools analysis --- 17.6.2022
########################################################################################
df_filt_crc = df_hg19_filtered_crc.sort_values(by=["Chr_Start"])
new_col = ["chromosome", "Chr_Start", "Chr_End", "symbol", "Exon", "Strand"]
df_filt_crc = df_filt_crc[new_col]
df_filt_crc["chromosome"] = df_filt_crc['chromosome'].apply(lambda x: f"chr{x}")

index_column = df_filt_crc.index
df_filt_crc=df_filt_crc.reset_index()
for i in range(len(df_filt_crc)-1):
    if not df_filt_crc["Strand"][i] == "-":
        continue
    else:
        df_filt_crc["Chr_Start"][i], df_filt_crc["Chr_End"][i] = df_filt_crc["Chr_End"][i], df_filt_crc["Chr_Start"][i]   
del df_filt_crc["index"]

###### ------------------------------------------------------------------- ######  
###### ------ 02.07.2022 - Additional overlap ratio in bedtool_output ---- ######
###### ------------------------------------------------------------------- ###### 

def str_in_as_merge_psi_extended(bedtool_output, df_psi):
    """df_psi: Is the downloaded PSI values for all patients in the TCGASpliceSeq database
    --> here for cancer patients COAD (coad_psi) and READ (read_psi) summarized as CRC
    bedtool_output: is the output from intersect of tral panel and reference genomehg38 to hg19
    with following code in the command line tool
    # bedtools intersect -wao -a /mnt/c/Users/cmk_8/OneDrive/data/str_full_Hg38ToHg19.bed 
    # -b /mnt/c/Users/cmk_8/OneDrive/data/df_annotation_bed.bed > /mnt/c/Users/cmk_8/OneDrive/data/str_ref_overlap_Hg38toHg19.bed"""
    
    # map_psi= bedtools_output_preprocess(bedtool_output, df_psi)
        
    del bedtool_output[5] # Remove placeholder column from tral_perf_panel "."
    del bedtool_output[7] # Remove chromosome column, as it appears double
    bedtool_output = bedtool_output.rename(columns={0: "chromosome", 1: "start", 2: "end", 3: "repeat_length", 4: "str", 6: "db_match",
                                      8: "start_as", 9: "end_as", 10: "symbol", 11: "exon", 12: "strand", 13: "bp_overlap"})

    # Here just exchange " str_in_as" with desired bedtools output
    map_psi = df_psi[df_psi.symbol.isin(bedtool_output.symbol)].reset_index()
    map_psi = map_psi.rename(columns={"exons":"exon"})
    
    # Split out the exon-values, divided by ":"
    # Make it a Series, Stack the values
    exon_series = map_psi['exon'].str.split(':').apply(pd.Series, 1).stack()
    # Get rid of the stack:
    # Drop the level to line up with the DataFrame
    exon_series.index = exon_series.index.droplevel(-1)
    
    # Make your series a dataframe 
    exon_df = pd.DataFrame(exon_series)
    
    # Delete the `exon` column from the DataFrame
    del map_psi['exon']
    
    # Join the psi DataFrame to `exon_df`
    map_psi_new = map_psi.join(exon_df)
    
    # Check out the new `df`and bring column "exon" to front
    map_psi_new= map_psi_new.rename(columns={0: "exon"})
    cols= map_psi_new.columns.tolist()
    cols= cols[-1:] + cols[:-1]
    df_map_psi = map_psi_new[cols]
    
    # Merging the psi & as/str dataframes
    df_map_psi["exon"] = df_map_psi.exon.apply(pd.to_numeric, errors='coerce')
    #as_in_str_merge = as_in_str.merge(df_map_psi, how="inner", on=["symbol", "exon"])
    str_in_as_merge_overlap = bedtool_output.merge(df_map_psi, how="inner", on=["symbol", "exon"])

    # save the file as .tsv for later analysis in R
    str_in_as_merge_overlap.to_csv("C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/str_in_as_merge_overlap.tsv", sep="\t", header=True)
    
    as_type_count= str_in_as_merge_overlap.groupby("splice_type")["splice_type"].count()
    print(as_type_count)
    map_psi_count= df_map_psi.groupby("splice_type")["splice_type"].count()
    print(map_psi_count)
    
    str_in_as_merge_overlap.groupby(["splice_type", "repeat_length"])["splice_type"].count().unstack(fill_value=0).plot.bar()
    return str_in_as_merge_overlap, df_map_psi

    
str_as_border_overlap = pd.read_csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/str_as_overlap_Hg38ToHg19.bed', header=None, sep="\t")     
crc_psi = pd.read_csv("C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/PSI_COAD_and_READ.tsv", sep="\t")
str_in_as_merge_crc_overlap, map_psi_crc_overlap = str_in_as_merge_psi_extended(str_as_border_overlap, crc_psi)
#str_in_as_merge_crc_overlap.to_csv("C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/str_in_as_merge_crc_overlap.tsv", sep="\t", header=True)
str_in_as_merge_crc_overlap.groupby(["splice_type", "repeat_length"])["repeat_length"].count().unstack(fill_value=0).plot.bar()








# ----------------------------------------------------- Started, but not further investigated --------------------------
# 02.07.2022
str_in_as_crc_overlap = pd.read_csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/str_in_as_merge_crc_overlap.tsv',sep="\t")
str_in_as_crc_overlap.columns.patient.str.replace("_", "-")
str_in_as_crc_overlap["bp_overlap"] = str_in_as_crc_overlap["bp_overlap"]+1
crc_psi = pd.read_csv("C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/PSI_COAD_and_READ.tsv",  sep="\t")

# In this way, melt the clinical and "STR_PSI" dataset to map
# clinical information for AS event in each exon of each patient
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%matplotlib inline

# df1 = TCGA412[1:100]
# df2 = crc_psi.head()
# df2 = pd.melt(crc_psi, id_vars=crc_psi.columns[0:10], var_name= "patient", value_name = "psi")
# df2["patient"] = pd.DataFrame(df2.patient.str.replace("_", "-"))
# df_mss = df3[df3.merge(df4, on="patient")["MSI"] == "MSS"]
# df_msi = df3[df3.merge(TCGA412, on="patient")["MSI"] == "MSI"]


df3 = pd.melt(str_in_as_crc_overlap, id_vars=str_in_as_crc_overlap.columns[0:24], var_name= "patient", value_name = "psi")
del df3["Unnamed: 0"]
df3["patient"] = pd.DataFrame(df3.patient.str.replace("_", "-"))

########################################################################################
# Load MSI data (TCGA412)
########################################################################################
TCGA412 = pd.read_csv('C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/Semester 6/data/TCGA412.csv', sep=",")
#TCGA412.astype("object")
TCGA412= pd.DataFrame(TCGA412).rename(columns={"case_submitter_id": "patient"})

df4 = TCGA412

df_msi=df3.merge(df4, on="patient").query("MSI == 'MSI'")
df_mss=df3.merge(df4, on="patient").query("MSI == 'MSS'")

df_msi.patient.nunique()
df_mss.patient.nunique()

fig = plt.figure(figsize=(20, 10))
ax = (
    df_mss
        .groupby("splice_type")
        .size()
        .sort_values()
        .plot(kind="bar")
)
ax.set(ylabel="# STR loci genotyped", yscale=("log"))
plt.show()

                                          
    
                                              
    
                                              
    
                                              
    
                                              
    
                                              
    
                                              