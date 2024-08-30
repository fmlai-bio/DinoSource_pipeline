# DinoSource: a comprehensive database of dinoflagellate genomic resources

![image](https://github.com/fmlai-bio/DinoSource_pipeline/blob/main/DinoSource.png)
## DinoSource
[DinoSource website](http://glab.hzau.edu.cn/dinosource) <br>
[DinoSource tutorial](http://glab.hzau.edu.cn/dinosource/tutorial) <br>
[Contact with us](http://glab.hzau.edu.cn/dinosource/aboutus/contact)

## What is DinoSource?
  Here, we completed the genome assembly of Amphidinium carterae using PcaBio HiFi long reads and Hi-C error correction, annotating 41,274 genes based on Iso-seq, with 43 chromosomes. We have aggregated 20 publicly available dinoflagellate genomes and integrated data from various omics fields (including epigenomics, 3D genomics, transcriptomics, proteomics, and 695 datasets), providing a systematic and comprehensive perspective for understanding the biological mechanisms of dinoflagellates. Additionally, it offers a suite of comprehensive analytical tools to explore functional genes related to evolutionary biology or biological control.<br>
# DinoSource-pipeline
### The DinoSource pipeline have following dependencies :
* [Python (3.8.8) and Anaconda](https://www.anaconda.com/) for creating a needed environment.
* [samtools (>1.9)](http://www.htslib.org/download/) for processing alignment results 
* [pairtools (v1.0.2)](https://github.com/open2c/pairtools) for processing sequencing data from a Hi-C experiment
* [3DMax (v1.0)](https://github.com/BDM-Lab/3DMax) for 3D structure reconstructing
* [cooler (v0.9.1)](https://github.com/open2c/cooler) for matrix aggregation and performing operation on .mcool file
* [BWA (v0.7.17)](https://github.com/lh3/bwa) for aligning sequencing data from different experiments to genome
* [MACS2 (v2.1.1)](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) for peak calling used to identify areas in the genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing or a OCR-sequencing(e.g. ATAC-Seq) experiment
* [fastp (v0.20.0)](https://github.com/OpenGene/fastp) providing fast all-in-one preprocessing (QC/adapters/trimming/filtering/splitting/merging...) for FastQ files.
* [BatMeth2 (v2.01)](https://github.com/GuoliangLi-HZAU/BatMeth2) for Bisulfite DNA Methylation Data Analysis with Indel-sensitive Mapping
* [HiGlass Browser (v1.11)](https://docs.higlass.io/)
* [WashU epigenome browser (v54.0.2)](http://epigenomegateway.wustl.edu/)
## Cite us

        

## Lab Homepage
[Welcome to Guoliang's Lab Website of Bioinformatics](http://glab.hzau.edu.cn/)
<br></br>
