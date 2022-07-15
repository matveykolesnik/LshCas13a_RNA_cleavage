WD="/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_expression_evaluation"
echo $WD
cd $WD

CPU=16

EcRefSeq="RefSeqs/Nontargeting_samples_refseqs.fasta"
LshRefSeq="RefSeqs/L_shahii_DSM_19757.fasta"

EcRefSeqIndexDir="RefSeqs/E_coli_nt_index"
LshRefSeqIndexDir="RefSeqs/L_shahii_index"

mkdir $EcRefSeqIndexDir
mkdir $LshRefSeqIndexDir

EcAlignmentsDir="Alignments/Shmakov_et_al/E_coli_RNASeq"
LshAlignmentsDir="Alignments/Shmakov_et_al/L_shahii_RNASeq"

mkdir $EcAlignmentsDir
mkdir $LshAlignmentsDir

EcRefSeqIndexPath=${EcRefSeqIndexDir}/"E_coli_nt_bowtie2_index"
LshRefSeqIndexPath=${LshRefSeqIndexDir}/"L_shahii_bowtie2_index"

bowtie2-build $EcRefSeq $EcRefSeqIndexPath
bowtie2-build $LshRefSeq $LshRefSeqIndexPath

EcTrimmedDataDir="Data/Shmakov_et_al_2015/E_coli_RNASeq/Trimmed"
LshTrimmedDataDir="Data/Shmakov_et_al_2015/L_shahii_RNASeq/Trimmed"

Ec_prefix=$(ls $EcTrimmedDataDir | sed 's/_paired_R[1,2].fastq.gz//g' | sed 's/_unpaired_R[1,2].fastq.gz//g' | uniq)
Lsh_prefix=$(ls $LshTrimmedDataDir | sed 's/_paired_R[1,2].fastq.gz//g' | sed 's/_unpaired_R[1,2].fastq.gz//g' | uniq)

bowtie2 -q -x $EcRefSeqIndexPath -1 ${EcTrimmedDataDir}/${Ec_prefix}_paired_R1.fastq.gz -2 ${EcTrimmedDataDir}/${Ec_prefix}_paired_R2.fastq.gz -p $CPU | pigz -9 > ${EcAlignmentsDir}/E_coli_data.sam.gz

bowtie2 -q -x $LshRefSeqIndexPath -1 ${LshTrimmedDataDir}/${Lsh_prefix}_paired_R1.fastq.gz -2 ${LshTrimmedDataDir}/${Lsh_prefix}_paired_R2.fastq.gz -p $CPU | pigz -9 > ${LshAlignmentsDir}/L_shahii_data.sam.gz
