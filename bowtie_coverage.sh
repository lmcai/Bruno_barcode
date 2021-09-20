#1. get reference sequence in 3_barcode_extraction
#convert multiple line fasta to single line fasta | grep the header

SPECIES=$1
READS=$2
GENE=$3
MIN_COVERAGE_DEPTH=1


awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' 3_barcode_extraction/${GENE}.${READS}.1.fas | grep -A 1 "${SPECIES}$" >5_base_coverage/${GENE}.${READS}.${SPECIES}.fas

#2. Reads mapping with Bowtie2
#module load bowtie2/2.3.2-fasrc02
#module load samtools/1.10-fasrc01

bowtie2-build 5_base_coverage/${GENE}.${READS}.${SPECIES}.fas 5_base_coverage/${GENE}.${READS}.${SPECIES}
bowtie2 -x 5_base_coverage/${GENE}.${READS}.${SPECIES} --no-unal -1 0_FASTQ/${SPECIES}_dedup_R1.${READS}.1.fastq -2 0_FASTQ/${SPECIES}_dedup_R2.${READS}.1.fastq -S - | \
samtools view -bS - | \
samtools sort -m 5G > 5_base_coverage/${GENE}.${READS}.${SPECIES}.sorted.bam
samtools index 5_base_coverage/${GENE}.${READS}.${SPECIES}.sorted.bam 

#3. calculate mean base coverage
samtools depth -a --reference 5_base_coverage/${GENE}.${READS}.${SPECIES}.fas 5_base_coverage/${GENE}.${READS}.${SPECIES}.sorted.bam | awk '{print $3}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' >${SPECIES}.${READS}.${GENE}.cov

