#mkdir 01fastp_clean 02hisat2_results 03stringtie logs

vim run.sh
work_dir="/public/home/xxx/data/RNA-seq/Ac-N"
genome="/public/home/xxx/data/Acgenome/Ac.v1.genome.fa"
STAR_index="/public/home/xxx/data/Acgenome/star_index/Ac.v1.genome"
gtf="/public/home/xxx/data/Acgenome/Amphidinium_carterae.v2.representative.gtf"
nthread=12
fastq=$1

cd ${work_dir}
echo "trim_galore start"
trim_galore -q 25 --phred33 --fastqc --stringency 3 --length 20 -e 0.1 --paired --gzip -o ./01trim_galore  00raw_data/${fastq}_1.fastq.gz 00raw_data/${fastq}_2.fastq.gz

echo "STAR start"
STAR --runThreadN ${nthread} --genomeDir ${STAR_index} --readFilesIn 01trim_galore/${fastq}_1_val_1.fq.gz 01trim_galore/${fastq}_2_val_2.fq.gz --readFilesCommand zcat --sjdbGTFfile ${gtf} --sjdbOverhang 100 --outFileNamePrefix 02STAR/${fastq}. --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts

echo "samtools start"
samtools flagstat -@ ${nthread}  02STAR/${fastq}.Aligned.sortedByCoord.out.bam > 02STAR/${fastq}.Aligned.sortedByCoord.out.bam.stats
samtools sort 02STAR/${fastq}.Aligned.sortedByCoord.out.bam -@ ${nthread} -o 02STAR/${fastq}.sort.bam
samtools index -@ ${nthread} 02STAR/${fastq}.sort.bam

echo "bam2bw"
bamCoverage -b 02STAR/${fastq}.sort.bam -p 12 -o 02STAR/${fastq}.bam.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 1258288479



#qsub
for i in `cat file.txt`;do echo "bash /public/home/xxx/data/RNA-seq/Ac-N/run.sh ${i}" |qsub -l nodes=1:ppn=12 -N ${i} -d ./logs -j eo ;done



gtf="/public/home/fmlai/2023-11Dinosource/data/Acgenome/Amphidinium_carterae.v2.representative.gtf"
#featureCounts
nohup featureCounts -T 20 -p -t exon -g gene_id -a ${gtf} -o counts_gtf.txt 02hisat2_results/*sort.bam > logs/featureCounts_gtf.log 2>&1 &
#clean
cat counts.txt|sed '1d'|cut -f 1,7-|sed 's/.sort.bam//g'|sed 's/02hisat2_results\///g'|sed 's/gene://g' > 03featureCount/counts2.txt
