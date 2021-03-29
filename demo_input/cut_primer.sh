#!/bin/bash
### Jianshu Zhao (jianshu.zhao@gatech.edu). Primer remove bash script based on cutadapt
for F in *.fastq.gz; do
	gunzip $F
done
mkdir tmp
for F in *_R1.fastq; do
	R=${F%_*}_R2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	cutadapt -j $(nproc) -a ^GTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCC -A ^GGACTACNVGGGTWTCTAAT...TTACCGCGGCKGCTGRCAC -m 200 --discard-untrimmed -o ./tmp/${SAMPLE}_R1.fastq -p ./tmp/${SAMPLE}_R2.fastq ${F} ${R} > ./tmp/${SAMPLE}.cutadaptStats.txt 
done
cd tmp
for F in *_R1.fastq; do
        R=${F%_*}_R2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	pigz -p $(nproc) $F
	pigz -p $(nproc) $R
	mv ${SAMPLE}_R1.fastq.gz ../
	mv ${SAMPLE}_R2.fastq.gz ../
done
cd ../
for F in *.fastq; do
	rm $F
done
rm -rf ./tmp
