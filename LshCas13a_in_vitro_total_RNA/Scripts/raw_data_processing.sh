WD="LshCas13a_in_vitro_total_RNA/"
cd $WD
CPU=20

mkdir QC
mkdir QC/Trimmed
mkdir Data/Trimmed

Adapters="Reference_sequences/adapters.fa"

fastqc -t $CPU -o QC Data/*.fastq.gz

for i in $(ls Data/*fastq.gz | sed 's/_R[1,2]_001.fastq.gz//g' | uniq);
do
    echo $i;
    prefix=`basename $i`
    trimmomatic PE -threads $CPU -phred33 $i"_R1_001.fastq.gz" $i"_R2_001.fastq.gz" Data/Trimmed/${prefix}_paired_R1.fastq.gz Data/Trimmed/${prefix}_unpaired_R1.fastq.gz Data/Trimmed/${prefix}_paired_R2.fastq.gz Data/Trimmed/${prefix}_unpaired_R2.fastq.gz ILLUMINACLIP:$Adapters:2:30:7 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:15;
done

fastqc -t $CPU -o QC/Trimmed Data/Trimmed/*.fastq.gz
