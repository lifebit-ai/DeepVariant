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




<p align="center">
  <img src="https://github.com/cthulhu-tech/DeepVariant/blob/master/dpwf.jpg">
</p>

## INPUT PARAMETERS

### About preprocessing

DeepVariant, in order to run at its fastest, requires some indexed and compressed versions of both the reference genome and the BAM files. With DeepVariant in Nextflow, if you wish, you can only use as an input the fasta and the BAM file and let us do the work for you in a clean and standarized way (standard tools like [samtools](http://samtools.sourceforge.net/) are used for indexing and every step is run inside of  a Docker container).

This is how the input files that are needed as input 
```
NA12878_S1.chr20.10_10p1mb.bam   NA12878_S1.chr20.10_10p1mb.bam.bai	
ucsc.hg19.chr20.unittest.fasta   ucsc.hg19.chr20.unittest.fasta.fai 
ucsc.hg19.chr20.unittest.fasta.gz  ucsc.hg19.chr20.unittest.fasta.gz.fai   ucsc.hg19.chr20.unittest.fasta.gz.gzi

```
If you do not have all of them, these are the file you can give as input to the Nextflow pipeline, and the rest will be produced for you.
```
NA12878_S1.chr20.10_10p1mb.bam  
ucsc.hg19.chr20.unittest.fasta
```

### PARAMETERS DEFINITION 

- ### REFERENCE GENOME



------------------
In the following lines an overview on how DeepVariants parameters from DeepVariant are called in this Nextflow Pipeline.

- **params.input** folder which contains all the needed files as input ( BAM, BAI, FASTA etc)
    Here an exact example of which file should be contained: 
    
    * Example of content:
     NA12878_S1.chr20.10_10p1mb.bam			ucsc.hg19.chr20.unittest.fasta
     NA12878_S1.chr20.10_10p1mb.bam.bai		ucsc.hg19.chr20.unittest.fasta.fai
     test_nist.b37_chr20_100kbp_at_10mb.bed		ucsc.hg19.chr20.unittest.fasta.gz
     test_nist.b37_chr20_100kbp_at_10mb.vcf.gz	ucsc.hg19.chr20.unittest.fasta.gz.fai
     test_nist.b37_chr20_100kbp_at_10mb.vcf.gz.tbi	ucsc.hg19.chr20.unittest.fasta.gz.gzi


- **params.modelFolder** folder which contains all files for the prediction model 
    * Example of content:
    model.ckpt.data-00000-of-00001	model.ckpt.index		model.ckpt.meta
    
- **params.modelName** name of the model to be used contained in the modelFolder 
    * Example :
    model.ckpt
    
- **params.bam_name** name of the bam file from the input folder to be used.
    * Example : 
    params.bam_name="NA12878_S1.chr20.10_10p1mb.bam";
    
- **params.ref_name** name of the reference fasta genome from the input folder to be used. 
    * Example : 
    params.ref_name="ucsc.hg19.chr20.unittest.fasta";
    
- **params.regions** regions which need to be analyzed. 
    * Example: 
    params.regions="chr20:10,000,000-10,010,000";




    
    

    
    
