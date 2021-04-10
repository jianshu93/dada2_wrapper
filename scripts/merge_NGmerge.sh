for F in *_fastq.gz; do
	gunzip $F
done
mkdir tmp
for F in *_R1.fastq; do
	R=${F%_*}_R2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	NGmerge -1 ${F} -2 ${R} -o ./tmp/${SAMPLE}_R1.fastq -n $(nproc) -m 20
done

for F in *_R1.fastq; do
    R=${F%_*}_R2.fastq
	BASE=${F##*/}
	SAMPLE=${BASE%_*}
	pigz -p $(nproc) $F
	pigz -p $(nproc) $R
done