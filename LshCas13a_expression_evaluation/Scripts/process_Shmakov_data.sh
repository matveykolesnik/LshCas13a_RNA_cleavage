WD="/home/ledum/data/LshCas13a_RNA_cleavage/LshCas13a_expression_evaluation"
echo $WD
cd $WD
CPU=16

EcRawDataDir="Data/Shmakov_et_al_2015/E_coli_RNASeq/Raw"
LshRawDataDir="Data/Shmakov_et_al_2015/L_shahii_RNASeq/Raw"

EcTrimmedDataDir="Data/Shmakov_et_al_2015/E_coli_RNASeq/Trimmed"
LshTrimmedDataDir="Data/Shmakov_et_al_2015/L_shahii_RNASeq/Trimmed"

EcRawDataQCDir="Data/Shmakov_et_al_2015/E_coli_RNASeq/RawQC"
EcTrimmedDataQCDir="Data/Shmakov_et_al_2015/E_coli_RNASeq/TrimmedQC"
LshRawDataQCDir="Data/Shmakov_et_al_2015/L_shahii_RNASeq/RawQC"
LshTrimmedDataQCDir="Data/Shmakov_et_al_2015/L_shahii_RNASeq/TrimmedQC"

mkdir $EcTrimmedDataDir
mkdir $LshTrimmedDataDir

mkdir $EcRawDataQCDir
mkdir $EcTrimmedDataQCDir
mkdir $LshRawDataQCDir
mkdir $LshTrimmedDataQCDir

Adapters="RefSeqs/Illumina_adapters.fasta"

fastqc -t $CPU -o $EcRawDataQCDir ${EcRawDataDir}/*.fastq.gz
fastqc -t $CPU -o $LshRawDataQCDir ${LshRawDataDir}/*.fastq.gz

Ec_prefix=$(ls $EcRawDataDir | sed 's/_[1,2].fastq.gz//g' | uniq)
Lsh_prefix=$(ls $LshRawDataDir | sed 's/_[1,2].fastq.gz//g' | uniq)


trimmomatic PE -threads $CPU -phred33 ${EcRawDataDir}/${Ec_prefix}_1.fastq.gz ${EcRawDataDir}/${Ec_prefix}_2.fastq.gz ${EcTrimmedDataDir}/${Ec_prefix}_paired_R1.fastq.gz ${EcTrimmedDataDir}/${Ec_prefix}_unpaired_R1.fastq.gz ${EcTrimmedDataDir}/${Ec_prefix}_paired_R2.fastq.gz ${EcTrimmedDataDir}/${Ec_prefix}_unpaired_R2.fastq.gz ILLUMINACLIP:$Adapters:2:30:7 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:15

trimmomatic PE -threads $CPU -phred33 ${LshRawDataDir}/${Lsh_prefix}_1.fastq.gz ${LshRawDataDir}/${Lsh_prefix}_2.fastq.gz ${LshTrimmedDataDir}/${Lsh_prefix}_paired_R1.fastq.gz ${LshTrimmedDataDir}/${Lsh_prefix}_unpaired_R1.fastq.gz ${LshTrimmedDataDir}/${Lsh_prefix}_paired_R2.fastq.gz ${LshTrimmedDataDir}/${Lsh_prefix}_unpaired_R2.fastq.gz ILLUMINACLIP:$Adapters:2:30:7 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:15

fastqc -t $CPU -o $EcTrimmedDataQCDir ${EcTrimmedDataDir}/*.fastq.gz
fastqc -t $CPU -o $LshTrimmedDataQCDir ${LshTrimmedDataDir}/*fastq.gz
