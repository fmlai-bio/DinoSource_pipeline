#mkdir 01batMeth2 02dmtools logs
vim run.sh
work_dir="/public/home/xxxx/data/WGBS/merge"
genome="/public/home/xxxx/data/Acgenome/batmeth2_index/Ac.v1.genome.fa"
nthread=24
fastq=$1

cd ${work_dir}
BatMeth2 align --fastp fastp -1 00raw_data/${fastq}_1.clean.fq.gz -2 00raw_data/${fastq}_2.clean.fq.gz -g $genome -o 01batMeth2/$fastq -p $nthread

dmtools bam2dm -C -E -S --Cx --remove_dup --maxcoverage 1000 --cf C -g $genome -b 01batMeth2/${fastq}.sort.bam -o 02dmtools/${fastq}.dm -p 1

dmtools ebsrate -i 02dmtools/${fastq}.dm -o 02dmtools/${fastq}.ebsrate.txt --bsmode chg

samtools flagstat 01batMeth2/${fastq}.sort.bam > 01batMeth2/${fastq}.sort.bam.stats



#qsub
for i in `cat file.txt`;do echo "bash /public/home/xxx/data/WGBS/merge/run.sh ${i}" |qsub -l nodes=1:ppn=24 -N ${i} -d ./logs -j eo ;done
