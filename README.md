# DinoSource: a comprehensive database of dinoflagellate genomic resources

![image](https://github.com/fmlai-bio/DinoSource_pipeline/blob/main/DinoSource.png)
## DinoSource
[DinoSource website](http://glab.hzau.edu.cn/dinosource) <br>
[DinoSource tutorial](http://glab.hzau.edu.cn/dinosource/tutorial) <br>
[Contact with us](http://glab.hzau.edu.cn/dinosource/aboutus/contact)

## What is DinoSource?
  Dinoflagellates, pivotal as primary producers in marine ecosystems and infamous for causing harmful algal blooms (red tides) that release toxins to threaten human health, profoundly impacting aquatic environments. They are a group of unicellular eukaryotes characterized by several distinctive features, such as permanently condensed liquid crystalline chromosomes without nucleosomes, unusually high amounts of 5-hydroxymethyluridine in their genomic DNA, abundant N1-methyladenosine modification in their mRNA, some of the largest genome sizes up to 80 times that of the human genome, and highly reduced plastids, providing valuable insights for biological and evolutionary studies.

Here, we present DinoSource — a comprehensive database for dinoflagellate genomics resources. DinoSource provides 21 genome assemblies that encompass all 20 currently sequenced dinoflagellate species, primarily representing 10 genus: Amoebophrya, Amphidinium, Breviolum, Cladocopium, Durusdinium, Effrenium, Fugacium, Polarella, Prorocentrum, and Symbiodinium. In addition, DinoSource integrates 703 omics samples, which cover various data types, including nucleotide modifications, DNA methylation, 3D genomics, chromatin accessibility, transcriptomics, and translatomics The database also offers a range of bioinformatics tools and genome browsers for intuitive data exploration and analysis.

DinoSource offers powerful tools for researchers to deeply understand the genomic structure and function of dinoflagellates, advancing research and applications in this field. We believe DinoSource will become an essential resource for exploring the ecological roles and evolutionary mechanisms of eukaryotes.<br>
# DinoSource-pipeline
### The DinoSource pipeline have following dependencies :
* [Python (3.8.8) and Anaconda](https://www.anaconda.com/) for creating a needed environment.
* [TrimGalore(0.6.10)](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for fastq data clean.
* [STAR (2.7.9a)](https://github.com/alexdobin/STAR) for RNA reads alignment
* [Bowtie2 (2.5.1)](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for DNA reads aligment
* [FeatureCounts (2.0.3)](https://subread.sourceforge.net/featureCounts.html) RNAseq quantitative tool
* [samtools (>1.9)](http://www.htslib.org/download/) for processing alignment results 
* [pairtools (v1.0.2)](https://github.com/open2c/pairtools) for processing sequencing data from a Hi-C experiment
* [3DMax (v1.0)](https://github.com/BDM-Lab/3DMax) for 3D structure reconstructing
* [cooler (v0.9.1)](https://github.com/open2c/cooler) for matrix aggregation and performing operation on .mcool file
* [BWA (v0.7.17)](https://github.com/lh3/bwa) for aligning sequencing data from different experiments to genome
* [MACS2 (v2.1.1)](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) for peak calling used to identify areas in the genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing or a OCR-sequencing(e.g. ATAC-Seq) experiment
* [fastp (v0.20.0)](https://github.com/OpenGene/fastp) providing fast all-in-one preprocessing (QC/adapters/trimming/filtering/splitting/merging...) for FastQ files.
* [BatMeth2 (v2.01)](https://github.com/GuoliangLi-HZAU/BatMeth2) for Bisulfite DNA Methylation Data Analysis with Indel-sensitive Mapping
* [HiGlass Browser (v1.11)](https://docs.higlass.io/)
* [WashU epigenome browser (v54.0.4)](http://epigenomegateway.wustl.edu/)
## Cite us

        

## Lab Homepage
[Welcome to Guoliang's Lab Website of Bioinformatics](http://glab.hzau.edu.cn/)
<br></br>
