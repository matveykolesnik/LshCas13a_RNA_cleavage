WD="/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_expression_evaluation/"
echo $WD
cd $WD

CPU=16
IS=200

EcC3000TrimmedDataDir="../LshCas13a_C3000/Data/Trimmed"
Ecd10LVMTrimmedDataDir="../LshCas13a_d10LVM/Data/Trimmed"

EcC3000BAMDir="Alignments/E_coli_C3000_BAM"
Ecd10LVMBAMDir="Alignments/E_coli_d10LVM_BAM"

mkdir $EcC3000BAMDir
mkdir $Ecd10LVMBAMDir

EcRef="RefSeqs/Nontargeting_samples_refseqs.fasta"
EcRefIndex="RefSeqs/Nontargeting_samples_refseqs_index"

bowtie2-build $EcRef $EcRefIndex

for i in $(ls ${EcC3000TrimmedDataDir}/*.fastq.gz | sed 's/_paired_R[1,2].fastq.gz//g' | sed 's/_unpaired_R[1,2].fastq.gz//g' | uniq);
do
	echo $i;
	prefix=$(basename $i);
	bowtie2 -x $EcRefIndex -1 ${i}"_paired_R1.fastq.gz" -2 ${i}"_paired_R2.fastq.gz" -X $IS -p $CPU | samtools view -Sbu - | samtools sort -@ 16 - -o ${EcC3000BAMDir}/${prefix}_sorted.bam;
	samtools index ${EcC3000BAMDir}/${prefix}_sorted.bam;

done

for i in $(ls ${Ecd10LVMTrimmedDataDir}/*.fastq.gz | sed 's/_paired_R[1,2].fastq.gz//g' | sed 's/_unpaired_R[1,2].fastq.gz//g' | uniq);
do
	echo $i;
	prefix=$(basename $i);
	bowtie2 -x $EcRefIndex -1 ${i}"_paired_R1.fastq.gz" -2 ${i}"_paired_R2.fastq.gz" -X $IS -p $CPU | samtools view -Sbu - | samtools sort -@ 16 - -o ${Ecd10LVMBAMDir}/${prefix}_sorted.bam;
	samtools index ${Ecd10LVMBAMDir}/${prefix}_sorted.bam;

done
