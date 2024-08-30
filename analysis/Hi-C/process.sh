#mkdir 01fastp_clean 02bwa 03pairtools 04cooler logs
vim hicrun.sh
work_dir="/public/home/xxx/data/HIC"
genome_index="/public/home/xxx/data/Acgenome/bwa_index/Ac.v1.genome"
genome_fai="/public/home/xxx/data/Acgenome/Ac.v1.genome.fa.fai"
bin_size=1000
nthreads=12
fastq=$1

cd ${work_dir}
export PATH=/public/home/xxx/miniconda3/envs/pairtools0.3.0/bin:$PATH
#qc
fastp -w ${nthreads} -i 00raw_data/${fastq}_1.fq.gz -I 00raw_data/${fastq}_2.fq.gz -o 01fastp_clean/${fastq}_1.clean.fq.gz -O 01fastp_clean/${fastq}_2.clean.fq.gz -h 01fastp_clean/${fastq}.fastp.html -j 01fastp_clean/${fastq}.fastp.json
#alignment
echo "1.align start"
bwa mem -SP5M -t $nthreads $genome_index 01fastp_clean/${fastq}_1.clean.fq.gz 01fastp_clean/${fastq}_2.clean.fq.gz |samtools view -hbS > 02bwa/${fastq}.bam
echo "samtools start"
samtools flagstat -@ ${nthreads} 02bwa/${fastq}.bam > 02bwa/${fastq}.bam.stats
echo "1.align end"
#filtering
echo "2.pairparse start"
/public/home/xxx/miniconda3/envs/pairtools0.3.0/bin/pairtools parse -c ${genome_fai} 02bwa/${fastq}.bam -o 03pairtools/${fastq}.pairsam.gz --drop-sam
echo "2.pairparse end"
echo "3.pairsort start"
/public/home/xxx/miniconda3/envs/pairtools0.3.0/bin/pairtools sort --memory 24G 03pairtools/${fastq}.pairsam.gz -o 03pairtools/${fastq}.sort.pairsam.gz
echo "3.pairsort end"
echo "4.pairdedup start"
/public/home/xxx/miniconda3/envs/pairtools0.3.0/bin/pairtools dedup --mark-dups 03pairtools/${fastq}.sort.pairsam.gz -o 03pairtools/${fastq}.dedup.pairsam.gz --output-stats 03pairtools/${fastq}.dedup.stats --output-dups 03pairtools/${fastq}.dup.pairsam.gz
echo "4.pairdedup end"
echo "5.pairselect end"
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' 03pairtools/${fastq}.dedup.pairsam.gz -o 03pairtools/${fastq}.select.pairsam.gz
echo "5.pairselect end"
echo "6.pairindex start"
pairix -p pairs 03pairtools/${fastq}.select.pairsam.gz
echo "6.pairindex end"
#matrix aggregation and normalization
echo "7.2cool start"
#cooler cload pairix ${genome_fai}:${bin_size} 03pairtools/${fastq}.select.pairsam.gz 04cooler/${fastq}.${bin_size}.cool
#or
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 ${genome_fai}:${bin_size} 03pairtools/${fastq}.select.pairsam.gz 04cooler/${fastq}.${bin_size}.cool
echo "7.2cool end"
cooler zoomify -p 5 -r 1000,2000,5000N --balance 04cooler/${fastq}.${bin_size}.cool -o 04cooler/${fastq}.mcool
echo "8.pairtools stats start"
pairtools stats 03pairtools/${fastq}.pairsam.gz -o 03pairtools/${fastq}.pairsam.gz.stats
#pairtools stats 03pairtools/${fastq}.dup.pairsam.gz -o 03pairtools/${fastq}.dup.pairsam.gz.stats
pairtools stats 03pairtools/${fastq}.select.pairsam.gz -o 03pairtools/${fastq}.select.pairsam.gz.stats


# Batch submission task
for i in `cat file.txt`;do echo "bash /public/home/xxx/data/HIC/hicrun.sh ${i}"|qsub -l nodes=1:ppn=12 -d ./logs -N ${i}_hic1 -j eo ;done

#plot heatmap
hicPlotMatrix -m all_merge.keep43chr.nobalance.mcool::resolutions/1000000 -out test.png --vMax 1000 --dpi 800 --colorMap YlOrRd
hicPlotMatrix -m all_merge.keep43chr.nobalance.mcool::resolutions/50000 -out test-chr1.pdf --vMax 10 --dpi 800 --colorMap Reds --region Acart_chr01

