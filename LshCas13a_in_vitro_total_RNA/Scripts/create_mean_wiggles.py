#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 14:05:16 2020

@author: matvey
"""

import pandas as pd
import os

WD = "/home/matvey/data/LshCas13a_RNA_cleavage/LshCas13a_in_vitro_total_RNA"
OutPutDir = "Results/WIG_files/aveCPM_ends_counts"
os.chdir(WD)

os.makedirs(OutPutDir)

N5ETableCPMFile = "Results/Tables/Merged_ends_counts/N5E_T_vs_NT_CPM_normalized.tsv.gz"
N5ETableCPM = pd.read_csv(N5ETableCPMFile, sep="\t")

def create_wig_files(CountsTable, Column, OutPutDir, prefix):
    Seq_list = sorted(list(set(CountsTable.SeqID.tolist())))
    
    with open(os.path.join(OutPutDir, prefix+"_forward_aveCPM.wig"), "w") as hForwardWig:
        wig_header = f'track type=wiggle_0 name="{prefix} forward"\n'
        hForwardWig.write(wig_header)
        for refseq in Seq_list:
            hForwardWig.write(f"fixedStep chrom={refseq} start=1 step=1\n")
            wig_record = CountsTable[(CountsTable.SeqID == refseq) & (CountsTable.Strand == "+")][Column].tolist()
            for v in wig_record:
                hForwardWig.write(f"{v}\n")
    
    with open(os.path.join(OutPutDir, prefix+"_reverse_aveCPM.wig"), "w") as hReverseWig:
        wig_header = f'track type=wiggle_0 name="{prefix} reverse"\n'
        hReverseWig.write(wig_header)
        for refseq in Seq_list:
            hReverseWig.write(f"fixedStep chrom={refseq} start=1 step=1\n")
            wig_record = CountsTable[(CountsTable.SeqID == refseq) & (CountsTable.Strand == "+")][Column].tolist()
            for v in wig_record:
                hReverseWig.write(f"{v}\n")

create_wig_files(CountsTable=N5ETableCPM, Column="T_aveCPM", OutPutDir=OutPutDir, prefix="in_vitro_total_RNA_targeting")
create_wig_files(CountsTable=N5ETableCPM, Column="NT_aveCPM", OutPutDir=OutPutDir, prefix="in_vitro_total_RNA_nontargeting")
