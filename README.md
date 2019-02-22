# DeepVariant as a Nextflow pipeline

A Nextflow pipeline for running the [Google DeepVariant variant caller](https://github.com/google/deepvariant).

At [Lifebit](https://lifebit.ai/?utm_campaign=documentation&utm_source=github&utm_medium=web) we developed this pipeline to ease and reduce cost for variant calling analyses. You can test the pipeline through our Platform: [Deploit](https://deploit.lifebit.ai/app/home). This allows you to run Deepvariant over cloud in a matter of a couple of clicks: and for single users our service is completely free! 

![](http://g.recordit.co/iyMUgBgw5E.gif)

Alternatively if you want to run the pipeline on the command line checkout [nf-core/deepvariant](https://github.com/nf-core/deepvariant) for a newer improved version of the pipeline!

Read more about DeepVariant in Nextflow in our [Blog post](https://blog.lifebit.ai/post/deepvariant/?utm_campaign=documentation&utm_source=github&utm_medium=web)


## What is DeepVariant and why in Nextflow?

The Google Brain Team in December 2017 released a [Variant Caller](https://www.ebi.ac.uk/training/online/course/human-genetic-variation-i-introduction/variant-identification-and-analysis/what-variant) based on DeepLearning: DeepVariant.

In practice, DeepVariant first builds images based on the BAM file, then it uses a DeepLearning image recognition approach to obtain the variants and eventually it converts the output of the prediction in the standard VCF format. 

DeepVariant as a Nextflow pipeline provides several advantages to the users. It handles automatically, through **preprocessing steps**, the creation of some extra needed indexed and compressed files which are a necessary input for DeepVariant, and which should normally manually be produced by the users.
Variant Calling can be performed at the same time on **multiple BAM files** and thanks to the internal parallelization of Nextflow no resources are wasted.
Nextflow's support of Docker allows to produce the results in a computational reproducible and clean way by running every step inside of a **Docker container**.
Moreover, you can easily run DeepVariant as a Nextflow pipeline in the **cloud** through the Lifebit platform and let it do the hard work of configuring, scheduling and deploying for you.

For more detailed information about DeepVariant please refer to: 
https://github.com/google/deepvariant
https://research.googleblog.com/2017/12/deepvariant-highly-accurate-genomes.html


## Dependencies

[Nextflow](https://www.nextflow.io/)
[Docker](https://www.docker.com/)

## Test Run 

Locally or through the Deploit platform run: 

```
git clone https://github.com/lifebit-ai/DeepVariant
cd DeepVariant
nextflow run main.nf --test
```

In this way, the **prepared data** on repository (under test data) are used for running DeepVariant and you can find the produced VCF files in the folder "RESULTS-DeepVariant".

The input of the pipeline can be eventually changed as explained in the "Input parameters" section.

## Quick Start

A typical run on **whole genome data** looks like this: 
```
git clone https://github.com/lifebit-ai/DeepVariant
cd DeepVariant
nextflow run main.nf --hg19 --bam_folder "s3://deepvariant-data/test-bam/"
```
In this case variants are called on the two bam files contained in the lifebit-test-data/bam s3 bucket. The hg19 version of the reference genome is used.
Two vcf files are produced and can be found in the folder "RESULTS-DeepVariant"


A typical run on **whole exome data** looks like this: 
```
git clone https://github.com/lifebit-ai/DeepVariant
cd DeepVariant
nextflow run main.nf --exome --hg19 --bam_folder myBamFolder --bed myBedFile"
```


## More about the pipeline 

As shown in the following picture, the worklow both contains **preprocessing steps** ( light blue ones ) and proper **variant calling steps** ( darker blue ones ).

Some input files ar optional and if not given, they will be automatically created for the user during the preprocessing steps. If these are given, the preprocessing steps are skipped. For more information about preprocessing, please refer to the "INPUT PARAMETERS" section.

The worklow **accepts one reference genome and multiple BAM files as input**. The variant calling for the several input BAM files will be processed completely indipendently and will produce indipendent VCF result files. The advantage of this approach is that the variant calling of the different BAM files can be parallelized internally by Nextflow and take advantage of all the cores of the machine in order to get the results at the fastest.


<p align="center">
  <img src="https://github.com/lifebit-ai/DeepVariant/blob/master/pics/pic_workflow.jpg">
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
If you do not have all of them, these are the file you can give as input to the Nextflow pipeline, and the rest will be automatically  produced for you .
```
NA12878_S1.chr20.10_10p1b.bam  
ucsc.hg19.chr20.unittest.fasta
```

### Parameters definition 

- ### BAM FILES 

```
--bam_folder "/path/to/folder/where/bam/files/are"            REQUIRED
--getBai "true"                                               OPTIONAL  (default: "false")
```
In case only some specific files inside the BAM folder should be used as input, a file prefix can be defined by: 
```
--bam_file_prefix MYPREFIX
```

All the BAM files on which the variant calling should be performed should be all stored in the same folder. If you already have the index files (BAI) they should be stored in the same folder and called with the same prefix as the correspoding BAM file ( e.g. file.bam and file.bam.bai ). 

**! TIP** 
All the input files can be used in s3 buckets too and the s3://path/to/files/in/bucket can be used instead of a local path.


- ### REFERENCE GENOME

 By default the hg19 version of the reference genome is used. If you want to use it, you do not have to pass anything.
 
 If you do not want to use the deafult version, here is how it works: 
 
 Two standard version of the genome ( hg19 and GRCh38.p10 ) are prepared with all their compressed and indexed file in a   lifebit s3 bucket.
 They can be used by using one of the flags:
 
 ```
 --hg19 (default) 
 --h38
 --hs37d5
 --grch37primary
 ```
 
 For testing purposes we provide the chr20 of the hg19 version of the genome, accessible by: 
 ```
 --hg19chr20
 ```
 Alternatively, a user can use an own reference genome version, by using the following parameters:

  ```
  --fasta "/path/to/myGenome.fa"                OPTIONAL
  --fai   "/path/to/myGenome.fa.fai"            OPTIONAL
  --fastagz "/path/to/myGenome.fa.gz"           OPTIONAL
  --gzfai  "/path/to/myGenome.fa.gz.fai"         OPTIONAL
  --gzi  "/path/to/myGenome.fa"                  OPTIONAL
  ```
If the optional parameters are not passed, they will be automatically be produced for you and you will be able to find them in the "preprocessingOUTPUT" folder.

- ### Exome data and Bed file

If you are running on exome data you need to prodive the --exome flag so that the right verison of the model will be used.
Moreover, you can provide a bed file.
```
nextflow run main.nf --exome --hg19 --bam_folder myBamFolder --bed myBedFile
```


### Advanced parameters options

- ### CPUS 

The **make_example** process can be internally parallelized and it can be defined how many cpus should be assigned to this process.
By default all the cpus of the machine are used.

```
--j 2          OPTIONAL (default: all)
```
- ### MODEL 

The trained model which is used by the **call_variants** process can be changed.
The default one is the 0.6.0 Version for the whole genome. So if that is what you want to use too, nothing needs to be changed.
If you want to access the version 0.6.0 for the whole exome model, you need to use the --exome flag.
```
nextflow run main.nf --exome --hg19 --bam_folder myBamFolder --bed myBedFile
```

In case you want to use another version of the model you can change it by: 
```
--modelFolder "s3://deepvariant-test/models"
--modelName   "model.ckpt"
```
The modelName parameter describes the name of the model that should be used among the ones found in the folder defined by the parameter modelFolder.  The model folder must contain 3 files, the list of which looks like this:

```
model.ckpt.data-00000-of-00001
model.ckpt.index
model.ckpt.meta
```









    
    

    
    
