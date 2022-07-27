WD="/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_d10LVM"
echo $WD
cd $WD

RawDataDir="Data/Raw"
RawDataFilteredDir="Data/RawBBDuckCleaned"
RawDataContaminatedDir="Data/RawBBDukContaminated"
BBDukStatsDir="Data/BBDukStats"

mkdir $RawDataFilteredDir
mkdir $RawDataContaminatedDir
mkdir $BBDukStatsDir

S_cerevisiae_refseq="Reference_sequences/S_cerevisiae.fna"

for i in $(ls ${RawDataDir}/*fastq.gz | sed 's/_R[1,2]_001.fastq.gz//g' | uniq);
do
    echo $i;
    prefix=$(basename $i)

    bbduk.sh -Xmx25g in1=${RawDataDir}/${prefix}"_R1_001.fastq.gz" in2=${RawDataDir}/${prefix}"_R2_001.fastq.gz" out1=${RawDataFilteredDir}/${prefix}"_R1_001.fastq.gz" out2=${RawDataFilteredDir}/${prefix}"_R2_001.fastq.gz" outm1=${RawDataContaminatedDir}/${prefix}"_R1_001.fastq.gz" outm2=${RawDataContaminatedDir}/${prefix}"_R2_001.fastq.gz" ref=$S_cerevisiae_refseq k=31 hdist=1 stats=${BBDukStatsDir}/${prefix}"_stats.txt";
done
