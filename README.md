# DeepVariant as a Nextflow pipeline

It permits to run DeepVariant in its whole functionalities and it eases preprocessing steps

## What is DeepVariant and why in Nextflow?

The Google Brain Team in December 2017 released a [Variant Caller](https://www.ebi.ac.uk/training/online/course/human-genetic-variation-i-introduction/variant-identification-and-analysis/what-variant) based on DeepLearning: DeepVariant.

In practice, DeepVariant first builds images based on the BAM file, then it uses a DeepLearning image recognition approach to obtain the variants and eventually it converts the output of the prediction in the standard VCF format. 


WHY NEXTFLOW?
Preprocessing embedded
Parallelization 
Easy control 
Containerization 


For more detailed information about DeepVariant please refer to: 
https://github.com/google/deepvariant
https://research.googleblog.com/2017/12/deepvariant-highly-accurate-genomes.html


## Quick run

On a aws machine run: 

```
git clone https://github.com/cthulhu-tech/DeepVariant
cd DeepVariant
nextflow run main.nf
```

In this way, the **prepared data** on s3 bucket are used for running DeepVariant and you can find the produced VCF files in the folder "RESULTS-DeepVariant".
To run this step can be useful for a user to **see how it looks like to run the completely funcitonal version of the pipeline**.
The input of the pipeline can be eventually changed as explained in the "Input parameters" section.

## The workflow 

As shown in the following picture, the worklow both contains **preprocessing steps** ( light blue ones ) and proper **variant calling steps** ( darker blue ones ).

Some input files ar optional and if not given, they will be automatically created for the user during the preprocessing steps. If these are given, the preprocessing steps are skipped. For more information about preprocessing, please refer to the "INPUT PARAMETERS" section.

The worklow **accepts one reference genome and multiple BAM files as input**. The variant calling for the several input BAM files will be processed completely indipendently and will produce indipendent VCF result files. The advantage of this approach is that the variant calling of the different BAM files can be parallelized internally by Nextflow and take advantage of all the cores of the machine in order to get the results at the fastest.


<p align="center">
  <img src="https://github.com/cthulhu-tech/DeepVariant/blob/master/dpwf.jpg">
</p>

## INPUT PARAMETERS

### About preprocessing

DeepVariant, in order to run at its fastest, requires some indexed and compressed versions of both the reference genome and the BAM files. With DeepVariant in Nextflow, if you wish, you can only use as an input the fasta and the BAM file and let us do the work for you in a clean and standarized way (standard tools like [samtools](http://samtools.sourceforge.net/) are used for indexing and every step is run inside of  a Docker container).

This is how the list of the needed input files looks like. If these are passed all as input parameters, the preprocessing steps will be skipped. 
```
NA12878_S1.chr20.10_10p1mb.bam   NA12878_S1.chr20.10_10p1mb.bam.bai	
ucsc.hg19.chr20.unittest.fasta   ucsc.hg19.chr20.unittest.fasta.fai 
ucsc.hg19.chr20.unittest.fasta.gz  ucsc.hg19.chr20.unittest.fasta.gz.fai   ucsc.hg19.chr20.unittest.fasta.gz.gzi

```
If you do not have all of them, these are the file you can give as input to the Nextflow pipeline, and the rest will be automatically  produced for you.
```
NA12878_S1.chr20.10_10p1b.bam  
ucsc.hg19.chr20.unittest.fasta
```

### PARAMETERS DEFINITION 

- ### REFERENCE GENOME

  ```
  --fasta "/path/to/myGenome.fa"                REQUIRED
  --fai   "/path/to/myGenome.fa.fai"            OPTIONAL
  --fastagz "/path/to/myGenome.fa.gz"           OPTIONAL
  -gzfai  "/path/to/myGenome.fa.gz.fai"         OPTIONAL
  -gzi  "/path/to/myGenome.fa"                  OPTIONAL
  ```


- ### BAM FILES 

```
--bam_folder "/path/to/folder/where/bam/files/are"            REQUIRED
--getBai "true"                                               OPTIONAL  (default: "false")
```

All the BAM files on which the variant calling should be performed should be all stored in the same folder. If you already have the index files (BAI) they should be stored in the same folder and called with the same prefix as the correspoding BAM file ( e.g. file.bam and file.bam.bai ). 


- ### REGIONS

```
--regions chr20                            all of chromosome 20
--regions chr20:10,000,000-11,000,000      10-11mb of chr20
--regions "chr20 chr21"                    chromosomes 20 and 21
```

You can refer to (here)[https://github.com/google/deepvariant/blob/r0.5/docs/deepvariant-details.md] to get an exact overview of the parameter regions, which is used exactly used as it is in the standard version of DeepVariant.



- ### SHARDS 

The **make_example** process can be internally parallelized and it can be defined how many cores should be assigned to this process.
By default 2 cores are used.

```
--n_shards 2          OPTIONAL (default: 2)
```
- ### MODEL 

The trained model which is used by the **call_variants** process can be changed.
The default one is 









    
    

    
    
