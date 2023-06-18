WD="/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection"
cd $WD

BAMFilesDir="Alignments/BAM"
BEDFilesDir="Results/BEDFiles"
bigwigFilesDir="Results/bigwigFilesDir"
mkdir $BEDFilesDir
mkdir $bigwigFilesDir

#for i in ${BAMFilesDir}/*bam;
#do
#	echo $i;
#	prefix=$(basename $i .bam);
#	bamCoverage --bam $i --outFileName ${BEDFilesDir}/${prefix}_forward.bed --outFileFormat bedgraph --filterRNAstrand reverse --binSize 1 --numberOfProcessors 16 --normalizeUsing CPM;
#	bamCoverage --bam $i --outFileName ${BEDFilesDir}/${prefix}_reverse.bed --outFileFormat bedgraph --filterRNAstrand forward --binSize 1 --numberOfProcessors 16 --normalizeUsing CPM;
#done

for i in ${BAMFilesDir}/*bam;
do
	echo $i;
	prefix=$(basename $i .bam);
	bamCoverage --bam $i --outFileName ${bigwigFilesDir}/${prefix}_forward.bw --outFileFormat bigwig --filterRNAstrand reverse --binSize 1 --numberOfProcessors 16 --normalizeUsing CPM;
	bamCoverage --bam $i --outFileName ${bigwigFilesDir}/${prefix}_reverse.bw --outFileFormat bigwig --filterRNAstrand forward --binSize 1 --numberOfProcessors 16 --normalizeUsing CPM;
done
