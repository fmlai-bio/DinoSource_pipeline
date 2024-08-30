#mkdir 01trim_galore 02bowtie2 03picard 04alignmentSieve 05macs2 logs
vim run.sh
work_dir="/public/home/xxx/data/ATAC-seq/test_script"
bowtie2_index="/public/home/xxx/data/Acgenome/bowtie2_index/Ac.v1.genome"
fastq=$1
nthreads=10

cd $work_dir
echo "01、trim_galore"
trim_galore --paired --phred33 -q 25 -e 0.1 --stringency 3 --fastqc 00raw_data/${fastq}_1.fastq.gz 00raw_data/${fastq}_2.fastq.gz  -o ./01trim_galore 

echo "02、bowtie2"
#filter multi-align and Q < 30
bowtie2 --very-sensitive-local -X 2000 --no-unal -x ${bowtie2_index} -1 01trim_galore/${fastq}_1_val_1.fq.gz  -2 01trim_galore/${fastq}_2_val_2.fq.gz --un-gz 02bowtie2/${fastq}.unmapped.fq.gz -p ${nthreads}|samtools view -hbS  -q 30 > 02bowtie2/${fastq}.bam 
samtools sort -@ ${nthreads}-O bam -o 02bowtie2/${fastq}.sort.bam 02bowtie2/${fastq}.bam

echo "03、picard"
java -jar ~/software/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=02bowtie2/${fastq}.sort.bam O=03picard/${fastq}_picard.rmdup.bam M=03picard/${fastq}.picard.log 
samtools index -@ ${nthreads} 03picard/${fastq}_picard.rmdup.bam

echo "04、alignmentSieve"
alignmentSieve -p ${nthreads} --ATACshift --bam  03picard/${fastq}_picard.rmdup.bam -o 04alignmentSieve/${fastq}.tmp.bam
samtools sort -@ ${nthreads} -O bam -o 04alignmentSieve/${fastq}.shifted.bam 04alignmentSieve/${fastq}.tmp.bam
samtools index -@ ${nthreads} 04alignmentSieve/${fastq}.shifted.bam
rm 04alignmentSieve/${fastq}.tmp.bam

echo "05、macs2 call peak"
macs2 callpeak -f BAM --nomodel -g 1258288479 -B --keep-dup all --cutoff-analysis -t 04alignmentSieve/${fastq}.shifted.bam -n 05macs2/${fastq} -p 0.05

echo "06、bam2bw"
bamCoverage -b 04alignmentSieve/${fastq}.shifted.bam -p ${nthreads} -o 04alignmentSieve/${fastq}.shifted.bam.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 1258288479
bamCoverage -b 03picard/${fastq}_picard.rmdup.bam -p ${nthreads} -o 03picard/${fastq}_picard.rmdup.bam.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 1258288479


#qsub
for i in `cat file.txt`;do echo bash /public/home/xxx/data/ATAC-seq/test_script/run.sh $i|qsub -l nodes=1:ppn=12 -N ${i} -d ./logs -j eo;done


#nohup
for i in `cat file.txt`;do echo "nohup bash /public/home/xxx/data/ATAC-seq/test_script/run.sh ${i} > logs/${i}.nohup.out 2>&1 &" ;done > qsub.sh
bash qsub.sh
