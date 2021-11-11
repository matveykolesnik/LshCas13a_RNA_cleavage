#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 22:43:44 2020

@author: matvey
"""

import pandas as pd
import gzip, os, re
from functools import reduce

WD = "/home/matvey/data/LshCas13a_RNA_cleavage/LshCas13a_C3000/WT_TSS_calling_old_assembly/"

EndsCountsDir = "Results/Tables/Ends_counts/"
MergedEndsCountsDir = "Results/Tables/Merged_ends_counts"
design_file = "design.tsv"

os.chdir(WD)
os.makedirs(MergedEndsCountsDir, exist_ok=True)

def merge_tables(files_list, EndsCountsDir, output_file):
    df_list = list()
    samples_list = list()
    
    for f in files_list:
        sample_name = re.search("(.+).N\dE.tsv.gz", f).group(1)
        samples_list.append(sample_name)
        print("Reading %s..." % f)
        
        with gzip.open(os.path.join(EndsCountsDir, f), "rt") as htable:
            EndsCountsDF = pd.read_csv(htable, sep="\t")
            
        EndsCountsDF.rename(columns={"Counts" : sample_name}, inplace=True)
        df_list.append(EndsCountsDF)
    
    print("Merging data frames...")
    MergedEndsCountsDF = reduce(lambda left, right: pd.merge(left, right, how="inner", on = ["SeqID", "Pos", "Strand"]), df_list)
    
    ordered_column_list = ["SeqID", "Pos", "Strand"] + samples_list
    MergedEndsCountsDF = MergedEndsCountsDF.reindex(columns=ordered_column_list)
    
    print("Writing merged data frame to file...")
    with gzip.open(output_file, "wt") as houtput:
        MergedEndsCountsDF.to_csv(houtput, sep="\t", index=False)
        
design_table = pd.read_csv(design_file, sep="\t")

##Merge N5E counts tables in TAP-untreated
#merge_tables(files_list = [s+".N5E.tsv.gz" for s in design_table[design_table.TAP_treatment == "untreated"].Sample.tolist()],
#             EndsCountsDir = EndsCountsDir,
#             output_file = os.path.join(MergedEndsCountsDir, "N5E_TAP_untreated_T_vs_NT.tsv.gz"))

#Merge N5E counts tables in nontargeting TAP_treated and TAP_untreated samples
merge_tables(files_list = [s+".N5E.tsv.gz" for s in design_table[design_table.Exp_group == "nontargeting"].Sample.tolist()],
             EndsCountsDir = EndsCountsDir,
             output_file = os.path.join(MergedEndsCountsDir, "N5E_NT_TAP_treated_vs_TAP_untreated.tsv.gz"))

##Merge N5E counts tables in targeting TAP_treated and TAP_untreated samples
#merge_tables(files_list = [s+".N5E.tsv.gz" for s in design_table[design_table.Exp_group == "targeting"].Sample.tolist()],
#             EndsCountsDir = EndsCountsDir,
#             output_file = os.path.join(MergedEndsCountsDir, "N5E_T_TAP_treated_vs_TAP_untreated.tsv.gz"))
