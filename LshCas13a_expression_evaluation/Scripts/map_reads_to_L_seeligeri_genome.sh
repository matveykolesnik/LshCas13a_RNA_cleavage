WD="/home/matvey/data/LshCas13a_RNA_cleavage/LshCas13a_expression_evaluation"
cd $WD

ReadsDir="Data/Meeske_et_al_2019/L_seeligeri_RNASeq_Raw"
SAMAlignmentsDir="Alignments/L_seeligeri/SAM"

L_seeligeri_refseq="RefSeqs/L_seeligeri_SLCC3954.fasta"
L_seeligeri_index_dir="RefSeqs/L_seeligeri_SLCC3954_index"
L_seeligeri_index=${L_seeligeri_index_dir}/index

mkdir $L_seeligeri_index_dir


bowtie2-build $L_seeligeri_refseq $L_seeligeri_index


for i in $(ls $ReadsDir | sed 's/_R[1,2].fastq.gz//g' | uniq);
do
	echo $i;
	bowtie2 -x $L_seeligeri_index -1 ${ReadsDir}/${i}_R1.fastq.gz -2 ${ReadsDir}/${i}_R2.fastq.gz -q -p 4 | pigz -1 > ${SAMAlignmentsDir}/${i}.sam.gz;

done
