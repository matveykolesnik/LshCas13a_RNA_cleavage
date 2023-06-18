WD="~/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection"
cd $WD
CPU=16

gzSAMFilesDir="Alignments/SAM"
BAMFilesDir="Alignments/BAM"

mkdir $BAMFilesDir

for i in ${gzSAMFilesDir}/*.sam.gz;
do
	echo $i;
	prefix=$(basename $i .sam.gz);
	zcat $i | samtools view -Sbu - | samtools sort -@ $CPU - -o ${BAMFilesDir}/${prefix}"_sorted.bam";
	samtools index ${BAMFilesDir}/${prefix}"_sorted.bam";
done
