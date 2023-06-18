#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 15:43:18 2019

@author: matvey
"""
import gzip, os
import pandas as pd
import numpy as np
from multiprocessing import Pool

#%% Path section
WD = "/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_5_min"
SamFilesDir = "Alignments/SAM/"

ReadsCoordsDir = "Results/Tables/Reads_coords/"
EndsCountsDir = "Results/Tables/Ends_counts/"
EndsCoverageDir = "Results/WIG_files/Ends_coverage/"

os.chdir(WD)
os.makedirs(ReadsCoordsDir, exist_ok=True)
os.makedirs(EndsCountsDir, exist_ok=True)
os.makedirs(EndsCoverageDir, exist_ok=True)

#%% SAM binary flags
read_unmapped_flag = 0b100
reversed_flag = 0b10000

#Dictionary with the lengths of reference sequences
RefSeq_lens = {"NC_000913.3" : 4641652,
               "pC003_Sp1" : 10403,
               "pC008" : 5778,
               "pC003_RFP_spacer" : 10403}

#%%
def create_reads_coords_table(sam_file, output_table):
    houtput_table = gzip.open(output_table, "wt")
    with gzip.open(sam_file, "rt") as hsamfile:
        for line in hsamfile:
            if line[0].startswith("@"):
                continue
            fields = line.rstrip().rsplit(sep="\t")
            if int(fields[1]) & read_unmapped_flag == read_unmapped_flag:
                continue
            SeqID = fields[2]
            if int(fields[1]) & reversed_flag == reversed_flag:
                read_start = int(fields[3]) + len(fields[9]) - 1
                read_end = int(fields[3])
                read_strand = "-"
            else:
                read_start = int(fields[3])
                read_end = int(fields[3]) + len(fields[9]) - 1
                read_strand = "+"
            houtput_table.write("\t".join([SeqID, str(read_start), str(read_end), read_strand]) + "\n")
            
    houtput_table.close()    
#%%
def create_ends_counts_table(ReadsCoordsFile, N5E_path):
    counts_file_fields = ["SeqID", "Pos", "Counts", "Strand"]
    ReadsCoordsDF = pd.read_csv(ReadsCoordsFile, sep="\t", names=["SeqID", "Start", "End", "Strand"])
    Seq_list = sorted(list(set(ReadsCoordsDF.SeqID.tolist())))
    N5E_records = list()
    
    for refseq in Seq_list:
        seqlen = RefSeq_lens[refseq]
        forward_N5E = np.zeros(seqlen)
        reverse_N5E = np.zeros(seqlen)
        
        for i, row in ReadsCoordsDF[ReadsCoordsDF.SeqID == refseq].iterrows():
            try:
                if row["Strand"] == "+":
                    forward_N5E[row["Start"]-1] += 1
                else:
                    reverse_N5E[row["Start"]-1] += 1
            except IndexError:
                print("Read %i..%i is out of range" % (row["Start"], row["End"]))
                continue
        print(f"{refseq}:\n\tForward strand N5E counts: {np.sum(forward_N5E)},\n\tReverse strand N5E counts: {np.sum(reverse_N5E)}")
        
        N5E_records.extend(list(zip([refseq for i in range(seqlen*2)], #SeqID
                                    [i+1 for i in range(seqlen)]*2, #Position
                                    np.concatenate([forward_N5E, reverse_N5E]), #Counts
                                    ["+" for i in range(seqlen)]+["-" for i in range(seqlen)]))) #Strand
        
    N5E_DF = pd.DataFrame.from_records(N5E_records, columns=counts_file_fields)
    
    with gzip.open(N5E_path, "wt") as hN5E_table:
        N5E_DF.to_csv(hN5E_table, sep="\t", index=False)

#%%
def create_ends_coverage_wigs(EndCountsFile, EndsCoverageDir):
    prefix = os.path.basename(EndCountsFile).rsplit(".", 2)[0]    
    EndsCoordsDF = pd.read_csv(EndCountsFile, sep="\t")
    Seq_list = sorted(list(set(EndsCoordsDF.SeqID.tolist())))
    
    CPM_norm_factor = sum(EndsCoordsDF.Counts)/10**6
    print("CPM factor: " + str(CPM_norm_factor))
    EndsCoordsDF.Counts = EndsCoordsDF.Counts/CPM_norm_factor
    
    #build "+"-strand WIG file
    with open(os.path.join(EndsCoverageDir, prefix+"_forward_CPM_norm.wig"), "w") as hForwardWig:
        wig_header = 'track type=wiggle_0 name="%s forward"\n' % prefix
        hForwardWig.write(wig_header)
        
        for refseq in Seq_list:
            hForwardWig.write("fixedStep chrom=%s start=1 step=1\n" % refseq)
            wig_record = EndsCoordsDF[(EndsCoordsDF.SeqID == refseq) & (EndsCoordsDF.Strand == "+")].Counts.tolist()
            for v in wig_record:
                hForwardWig.write("%f\n" % v)
    
    with open(os.path.join(EndsCoverageDir, prefix+"_reverse_CPM_norm.wig"), "w") as hReverseWig:
        wig_header = 'track type=wiggle_0 name="%s reverse"\n' % prefix
        hReverseWig.write(wig_header)
        
        for refseq in Seq_list:
            hReverseWig.write("fixedStep chrom=%s start=1 step=1\n" % refseq)
            wig_record = EndsCoordsDF[(EndsCoordsDF.SeqID == refseq) & (EndsCoordsDF.Strand == "-")].Counts.tolist()
            for v in wig_record:
                hReverseWig.write("%f\n" % v)

#%% Create TAB files
SamFilesList = sorted([f for f in os.listdir(SamFilesDir) if f.endswith(".sam.gz")])
ArgsList = [(os.path.join(SamFilesDir, f), 
             os.path.join(ReadsCoordsDir, os.path.basename(f).rsplit(".")[0]+".TAB.tsv.gz")) 
            for f in SamFilesList]
if __name__ == '__main__':
    pool = Pool()
    pool.starmap(create_reads_coords_table, ArgsList)
    pool.close()
    pool.join()
#%% Create ends counts tables
ReadsCoordsList = sorted([f for f in os.listdir(ReadsCoordsDir) if f.endswith(".tsv.gz")])
ArgsList = [(os.path.join(ReadsCoordsDir, f),
            os.path.join(EndsCountsDir, os.path.basename(f).rsplit(".")[0]+".N5E.tsv.gz"))
            for f in ReadsCoordsList]
if __name__ == '__main__':
    pool = Pool()
    pool.starmap(create_ends_counts_table, ArgsList)
    pool.close()
    pool.join()
#%% Create ends counts coverage files
EndsCountsFilesList = sorted([f for f in os.listdir(EndsCountsDir) if f.endswith(".tsv.gz")])
ArgsList = [(os.path.join(EndsCountsDir, f),
             EndsCoverageDir)
            for f in EndsCountsFilesList]
if __name__ == '__main__':
    pool = Pool()
    pool.starmap(create_ends_coverage_wigs, ArgsList)
    pool.close()
    pool.join()