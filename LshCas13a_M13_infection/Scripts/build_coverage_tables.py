#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 03:02:42 2022

@author: ledum
"""

import gzip, os, re
import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import reduce

#%% Path section
WD = "/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection/"
ReadsCoordsDir = "Results/Tables/Reads_coords/"
CoverageTablesDir = "Results/Tables/CoverageTables"
WIGCoverageDir = "Results/WIG_files/Reads_coverage"

os.chdir(WD)
os.makedirs(CoverageTablesDir, exist_ok=True)
os.makedirs(WIGCoverageDir, exist_ok=True)
#%% Reference sequence lengths
RefSeq_lens = {"NC_000913.3" : 4641652,
                "pC002" : 10606,
                "pC003_MS2mat_spacer" : 10403,
                "M13+Sp1rev_C-PFS" : 7354}
#%% 
def create_coverage_files(ReadsCoordsFile, CoverageFile):
    counts_file_fields = ["SeqID", "Pos", "Counts", "Strand"]
    ReadsCoordsDF = pd.read_csv(ReadsCoordsFile, sep="\t", names=["SeqID", "Start", "End", "Strand"])
    Seq_list = sorted(list(set(ReadsCoordsDF.SeqID.tolist())))
    
    coverage_records = []
    
    for refseq in Seq_list:
        seqlen = RefSeq_lens[refseq]
        forward_cov = np.zeros(seqlen)
        reverse_cov = np.zeros(seqlen)
        
        for index, row in ReadsCoordsDF[ReadsCoordsDF.SeqID == refseq].iterrows():
            try:
                if row["Strand"] == "+":
                    for i in range(row["Start"]-1, row["End"]):
                        forward_cov[i] += 1
                else:
                    for i in range(row["End"]-1, row["Start"]):
                        reverse_cov[i] += 1
            except IndexError:
                print("Read %i..%i is out of range" % (row["Start"], row["End"]))
                continue
            
        coverage_records.extend(list(zip([refseq for i in range(seqlen*2)], #SeqID
                                    [i+1 for i in range(seqlen)]*2, #Position
                                    np.concatenate([forward_cov, reverse_cov]), #Counts
                                    ["+" for i in range(seqlen)]+["-" for i in range(seqlen)]))) #Strand
        
    CoverageDF = pd.DataFrame.from_records(coverage_records, columns=counts_file_fields)
    with gzip.open(CoverageFile, "wt") as hCoverageFile:
        CoverageDF.to_csv(hCoverageFile, sep="\t", index=False)
#%%
def create_ends_coverage_wigs(EndCountsFile, EndsCoverageDir):
    prefix = os.path.basename(EndCountsFile).rsplit(".", 2)[0]    
    EndsCoordsDF = pd.read_csv(EndCountsFile, sep="\t")
    Seq_list = sorted(list(set(EndsCoordsDF.SeqID.tolist())))
    
    #build "+"-strand WIG file
    with open(os.path.join(EndsCoverageDir, prefix+"_forward.wig"), "w") as hForwardWig:
        wig_header = 'track type=wiggle_0 name="%s forward"\n' % prefix
        hForwardWig.write(wig_header)
        
        for refseq in Seq_list:
            hForwardWig.write("fixedStep chrom=%s start=1 step=1\n" % refseq)
            wig_record = EndsCoordsDF[(EndsCoordsDF.SeqID == refseq) & (EndsCoordsDF.Strand == "+")].Counts.tolist()
            for v in wig_record:
                hForwardWig.write("%f\n" % v)
    
    with open(os.path.join(EndsCoverageDir, prefix+"_reverse.wig"), "w") as hReverseWig:
        wig_header = 'track type=wiggle_0 name="%s reverse"\n' % prefix
        hReverseWig.write(wig_header)
        
        for refseq in Seq_list:
            hReverseWig.write("fixedStep chrom=%s start=1 step=1\n" % refseq)
            wig_record = EndsCoordsDF[(EndsCoordsDF.SeqID == refseq) & (EndsCoordsDF.Strand == "-")].Counts.tolist()
            for v in wig_record:
                hReverseWig.write("%f\n" % v)
#%%
ReadsCoordsList = sorted([f for f in os.listdir(ReadsCoordsDir) if f.endswith(".tsv.gz")])
ArgsList = [(os.path.join(ReadsCoordsDir, f),
            os.path.join(CoverageTablesDir, os.path.basename(f).rsplit(".")[0]+"_coverage.tsv.gz"))
            for f in ReadsCoordsList]
if __name__ == '__main__':
    pool = Pool()
    pool.starmap(create_coverage_files, ArgsList)
    pool.close()
    pool.join()
#%%
CoverageTablesFilesList = sorted([f for f in os.listdir(CoverageTablesDir) if f.endswith(".tsv.gz")])
ArgsList = [(os.path.join(CoverageTablesDir, f),
             WIGCoverageDir)
            for f in CoverageTablesFilesList]
if __name__ == '__main__':
    pool = Pool()
    pool.starmap(create_ends_coverage_wigs, ArgsList)
    pool.close()
    pool.join()

#%%
def merge_tables(files_list, CoverageCountsDir, output_file):
    df_list = list()
    samples_list = list()
    
    for f in files_list:
        sample_name = re.search("(.+)_coverage.tsv.gz", f).group(1)
        samples_list.append(sample_name)
        print("Reading %s..." % f)
        
        with gzip.open(os.path.join(CoverageCountsDir, f), "rt") as htable:
            CoverageCountsDF = pd.read_csv(htable, sep="\t")
            
        CoverageCountsDF.rename(columns={"Counts" : sample_name}, inplace=True)
        df_list.append(CoverageCountsDF)
    
    print("Merging data frames...")
    MergedCoverageCountsDF = reduce(lambda left, right: pd.merge(left, right, how="inner", on = ["SeqID", "Pos", "Strand"]), df_list)
    
    ordered_column_list = ["SeqID", "Pos", "Strand"] + samples_list
    MergedCoverageCountsDF = MergedCoverageCountsDF.reindex(columns=ordered_column_list)
    
    print("Writing merged data frame to file...")
    with gzip.open(output_file, "wt") as houtput:
        MergedCoverageCountsDF.to_csv(houtput, sep="\t", index=False)
#%%
design_file = "design.tsv"
MergedCoverageCountsDir = "Results/Tables/MergedCoverageCounts"
os.makedirs(MergedCoverageCountsDir, exist_ok=True)

design_table = pd.read_csv(design_file, sep="\t")
design_table.sort_values(by="Exp_group", ascending=False, inplace=True)

#Merge N5E counts tables in TAP-untreated
merge_tables(files_list = [s+"_coverage.tsv.gz" for s in design_table.Sample.tolist()],
             CoverageCountsDir = CoverageTablesDir,
             output_file = os.path.join(MergedCoverageCountsDir, "Reads_coverage_merged.tsv.gz"))