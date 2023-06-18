WD="/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_M13_infection"
echo $WD
cd $WD
CPU=16

RawDataDir="Data/Raw"

RawDataQCDir="Data/RawQC"
TrimmedDataDir="Data/Trimmed"
TrimmedDataQCDir="Data/TrimmedQC"

mkdir $RawDataQCDir
mkdir $TrimmedDataDir
mkdir $TrimmedDataQCDir

Adapters="Reference_sequences/adapters.fa"

#fastqc -t $CPU -o $RawDataQCDir ${RawDataDir}/*.fastq.gz

for i in $(ls ${RawDataDir}/*fastq.gz | sed 's/_R[1,2]_001.fastq.gz//g' | uniq);
do
    prefix=`basename $i`
    trimmomatic SE -threads $CPU -phred33 ${RawDataDir}/${prefix}"_R1_001.fastq.gz" ${TrimmedDataDir}/${prefix}_trimmed_R1.fastq.gz ILLUMINACLIP:$Adapters:2:30:7 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:15;
done

fastqc -t $CPU -o ${TrimmedDataQCDir} ${TrimmedDataDir}/*.fastq.gz