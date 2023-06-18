WD="/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection"
cd $WD

#Count mapped reads for each sample
for f in Alignments/BAM/*bam; 
do
	printf "%s\t%s\n" `basename $f` `samtools view -F 4 -c $f` >> Results/Tables/reads_per_lib.tsv; 
done
