WD="/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection"
echo $WD
cd $WD

CPU=16

nt_refs_path="Reference_sequences/M13_exp_nontargeting.fasta"
t_refs_path="Reference_sequences/M13_exp_targeting.fasta"

mkdir Reference_sequences/Nontargeting_samples_refseqs_index
mkdir Reference_sequences/Targeting_samples_refseqs_index

nt_index_path="Reference_sequences/Nontargeting_samples_refseqs_index/nt_index"
t_index_path="Reference_sequences/Targeting_samples_refseqs_index/t_index"

bowtie2-build $nt_refs_path $nt_index_path
bowtie2-build $t_refs_path $t_index_path

mkdir Alignments
mkdir Alignments/SAM

for i in $(cat design.tsv | grep -w "nontargeting" | cut -f 1);
do
    echo $i;
    bowtie2 -q -x $nt_index_path -U Data/Trimmed/${i}_trimmed_R1.fastq.gz -p $CPU | gzip > Alignments/SAM/`basename $i`.sam.gz;
done

for i in $(cat design.tsv | grep -w "targeting" | cut -f 1);
do
    echo $i;
    bowtie2 -q -x $t_index_path -U Data/Trimmed/${i}_trimmed_R1.fastq.gz -p $CPU | gzip > Alignments/SAM/`basename $i`.sam.gz;
done
