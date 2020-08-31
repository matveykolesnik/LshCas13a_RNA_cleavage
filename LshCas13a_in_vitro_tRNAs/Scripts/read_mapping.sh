WD="LshCas13a_in_vitro_tRNAs"
echo $WD
cd $WD

CPU=20
IS=200

refseq_path="Reference_sequences/E_coli_MS2_small_RNAs.fasta"

mkdir Reference_sequences/E_coli_MS2_small_RNAs_index
index_path="Reference_sequences/E_coli_MS2_small_RNAs_index/E_coli_MS2_small_RNAs_index"

bowtie2-build $refseq_path $index_path

mkdir Alignments
mkdir Alignments/SAM

for i in $(ls Data/Trimmed/*fastq.gz | sed 's/_paired_R[1,2].fastq.gz//g' | sed 's/_unpaired_R[1,2].fastq.gz//g' | uniq);
do
    echo $i;
    bowtie2 -q -x $index_path -1 ${i}_paired_R1.fastq.gz -2 ${i}_paired_R2.fastq.gz -X $IS -p $CPU | gzip > Alignments/SAM/`basename $i`.sam.gz;
done
