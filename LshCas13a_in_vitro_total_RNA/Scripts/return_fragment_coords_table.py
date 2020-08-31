#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 15:43:18 2019

@author: matvey
"""
import gzip, os
import pandas as pd

# Path section
WD = "LshCas13a_in_vitro_total_RNA/"
SamFilesDir = "Alignments/SAM/"

TabFilesDir = "Results/Tables/Read_pairs_TABs/"
FragCoordsDir = "Results/Tables/Fragments_coords/"
EndsCountsDir = "Results/Tables/Ends_counts/"
EndsCoverageDir = "Results/WIG_files/Ends_coverage/"

os.chdir(WD)
os.makedirs(TabFilesDir, exist_ok=True)
os.makedirs(FragCoordsDir, exist_ok=True)
os.makedirs(EndsCountsDir, exist_ok=True)
os.makedirs(EndsCoverageDir, exist_ok=True)

# SAM binary flags
proper_pair_flag = 0b10
reversed_flag = 0b10000

#Dictionary with the lengths of reference sequences
RefSeq_lens = {"NC_000913.3" : 4641652,
               "NC_001417.2" : 3569,
               "crRNA_precursor" : 58,
               "target_RNA" : 122}

#looping through SAM file records, selecting properly aligned read pairs and writing the coordinates and flags to separate file
def create_read_pairs_TAB(sam_file, tab_file):
    htabfile = gzip.open(tab_file, "wt")
    with gzip.open(sam_file, "rt") as hsamfile:
        for line in hsamfile:
            if line[0] == "@":
                continue #skip SAM file header
            fields = line.rsplit(sep="\t")
            if int(fields[1]) & proper_pair_flag == proper_pair_flag: #select only reads that form proper pairs
                htabfile.write("\t".join([fields[1], fields[2], fields[3], fields[8]])+"\n")
    htabfile.close()

#loop through reads pairs and return list of 0-based fragments coords with the direction
def create_fragments_coords_table(tab_file, output_table):
    houtput_table = gzip.open(output_table, "wt")
    with gzip.open(tab_file, "rt") as htabfile:
        for reads_pair in zip(*[iter(htabfile)]*2):
            R1 = reads_pair[0].rstrip().rsplit("\t") #first read in the pair
            R2 = reads_pair[1].rstrip().rsplit("\t") #second read in the pair
            
            if R1[1] == R2[1]: #both reads in the pair should be aligned to the same sequence
                SeqID = R1[1]
            else:
                continue
            
            if int(R1[0]) & reversed_flag == reversed_flag: #if the first read is alinged in reverse orientation => fragment is aligned in reverse orientation
                fragment_start = int(R2[2])
                fragment_end = fragment_start+int(R2[3])-1
                fragment_strand = "-"
            else:
                fragment_start = int(R1[2])
                fragment_end = fragment_start+int(R1[3])-1
                fragment_strand = "+"
            houtput_table.write("\t".join([SeqID, str(fragment_start), str(fragment_end), fragment_strand])+"\n")
    houtput_table.close()
    
def create_ends_counts_table(FragCoordsFile, N5E_path, N3E_path):
    counts_file_fields = ["SeqID", "Pos", "Counts", "Strand"]
    FragCoordsDF = pd.read_csv(FragCoordsFile, sep="\t", names=["SeqID", "Start", "End", "Strand"])
    Seq_list = sorted(list(set(FragCoordsDF.SeqID.tolist())))
    N5E_records = list()
    N3E_records = list()
    
    for refseq in Seq_list:
        seqlen = RefSeq_lens[refseq]
        forward_N5E = [0 for i in range(seqlen)]
        forward_N3E = [0 for i in range(seqlen)]
        reverse_N5E = [0 for i in range(seqlen)]
        reverse_N3E = [0 for i in range(seqlen)]
        print("Counting ends of fragments mapped to %s" % refseq)
        for i, row in FragCoordsDF[FragCoordsDF.SeqID == refseq].iterrows():
            try:
                if row["Strand"] == "+":
                    forward_N5E[row["Start"]-1] += 1
                    forward_N3E[row["End"]-1] += 1
                else:
                    reverse_N5E[row["End"]-1] += 1
                    reverse_N3E[row["Start"]-1] += 1
            except IndexError:
                print("Fragment %i..%i is out of range" % (row["Start"], row["End"]))
                continue
        
        N5E_records.extend(list(zip([refseq for i in range(seqlen*2)], #SeqID
                                    [i+1 for i in range(seqlen)]*2, #Position
                                    forward_N5E+reverse_N5E, #Counts
                                    ["+" for i in range(seqlen)]+["-" for i in range(seqlen)]))) #Strand
        
        N3E_records.extend(list(zip([refseq for i in range(seqlen*2)],
                                    [i+1 for i in range(seqlen)]*2,
                                    forward_N3E+reverse_N3E,
                                    ["+" for i in range(seqlen)]+["-" for i in range(seqlen)])))
        
    N5E_DF = pd.DataFrame.from_records(N5E_records, columns=counts_file_fields)
    N3E_DF = pd.DataFrame.from_records(N3E_records, columns=counts_file_fields)
    
    with gzip.open(N5E_path, "wt") as hN5E_table:
        N5E_DF.to_csv(hN5E_table, sep="\t", index=False)
        
    with gzip.open(N3E_path, "wt") as hN3E_table:
        N3E_DF.to_csv(hN3E_table, sep="\t", index=False)

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

SamFilesList = sorted([f for f in os.listdir(SamFilesDir) if f.endswith(".sam.gz")])
for f in SamFilesList:
    print("Processing SAM file %s..." % f)
    tab_file_path = os.path.join(TabFilesDir, os.path.basename(f).rsplit(".")[0]+".TAB.tsv.gz")
    
    print("Creating TAB file...")
    create_read_pairs_TAB(sam_file=os.path.join(SamFilesDir, f), 
                          tab_file=tab_file_path)


TabFilesList = sorted([f for f in os.listdir(TabFilesDir) if f.endswith(".tsv.gz")])
for f in TabFilesList:
    print("Processing TAB file %s..." % f)
    frag_coords_file_path = os.path.join(FragCoordsDir, os.path.basename(f).rsplit(".")[0]+".FragmentCoords.tsv.gz")
    print("Creating fragments coords table..")
    
    create_fragments_coords_table(tab_file=os.path.join(TabFilesDir, f),
                                  output_table=frag_coords_file_path)

FragCoordsList = sorted([f for f in os.listdir(FragCoordsDir) if f.endswith(".tsv.gz")])
for f in FragCoordsList:
    print("Processing fragments coords file %s..." % f)
    N5E_file = os.path.join(EndsCountsDir, os.path.basename(f).rsplit(".")[0]+".N5E.tsv.gz")
    N3E_file = os.path.join(EndsCountsDir, os.path.basename(f).rsplit(".")[0]+".N3E.tsv.gz")
    
    create_ends_counts_table(FragCoordsFile=os.path.join(FragCoordsDir, f),
                             N5E_path=N5E_file,
                             N3E_path=N3E_file)
    
EndsCountsFilesList = sorted([f for f in os.listdir(EndsCountsDir) if f.endswith(".tsv.gz")])
for f in EndsCountsFilesList:
    print("Processing ends counts file %s..." % f)
    create_ends_coverage_wigs(EndCountsFile=os.path.join(EndsCountsDir, f), EndsCoverageDir=EndsCoverageDir)
