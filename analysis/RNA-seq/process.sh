#mkdir 01fastp_clean 02hisat2_results 03stringtie logs

vim run.sh
work_dir="/public/home/xxx/data/RNA-seq/Ac-N"
genome="/public/home/xxx/data/Acgenome/Ac.v1.genome.fa"
hisat2_index="/public/home/xxx/data/Acgenome/hisat2_index/Ac.v1.genome"
gtf="/public/home/xxx/data/Acgenome/Amphidinium_carterae.v2.representative.gtf"
nthread=12
fastq=$1

cd ${work_dir}
echo "fastp start"
fastp -w ${nthread} -i 00raw_data/${fastq}_1.fq.gz -I 00raw_data/${fastq}_2.fq.gz -o 01fastp_clean/${fastq}_1.fastq_clean.gz -O 01fastp_clean/${fastq}_2.fastq_clean.gz -h 01fastp_clean/${fastq}.fastp.html -j 01fastp_clean/${fastq}.fastp.json 2>&1 | tee ${fastq}.fastp.log 

echo "hisat2 start"
hisat2 -p ${nthread} -x $hisat2_index -1 $work_dir/01fastp_clean/${fastq}_1.fastq_clean.gz -2 $work_dir/01fastp_clean/${fastq}_2.fastq_clean.gz |samtools view -hbS > $work_dir/02hisat2_results/${fastq}.bam

echo "samtools start"
samtools flagstat -@ ${nthread} $work_dir/02hisat2_results/${fastq}.bam > ./02hisat2_results/${fastq}.bam.stats
samtools sort ./02hisat2_results/${fastq}.bam -@ ${nthread} -o ./02hisat2_results/${fastq}.sort.bam
samtools index -@ ${nthread} ./02hisat2_results/${fastq}.sort.bam


echo "stringtie start"
stringtie 02hisat2_results/${fastq}.sort.bam -e -p ${nthread} -o ./03stringtie/${fastq}_transcript.gtf -G ${gtf}



#qsub
for i in `cat file.txt`;do echo "bash /public/home/xxx/data/RNA-seq/Ac-N/run.sh ${i}" |qsub -l nodes=1:ppn=12 -N ${i} -d ./logs -j eo ;done



gtf="/public/home/fmlai/2023-11Dinosource/data/Acgenome/Amphidinium_carterae.v2.representative.gtf"
#featureCounts
nohup featureCounts -T 20 -p -t exon -g gene_id -a ${gtf} -o counts_gtf.txt 02hisat2_results/*sort.bam > logs/featureCounts_gtf.log 2>&1 &
#clean
cat counts.txt|sed '1d'|cut -f 1,7-|sed 's/.sort.bam//g'|sed 's/02hisat2_results\///g'|sed 's/gene://g' > 03featureCount/counts2.txt
