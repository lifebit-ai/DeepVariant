# DeepVariant

Deep Variant as a Nextflow pipeline

For information about DeepVariant: 
https://github.com/google/deepvariant
https://research.googleblog.com/2017/12/deepvariant-highly-accurate-genomes.html


## Quick test: 

On a aws machine run: 

```
nextflow run main.nf
```

In this ways the prepared data on s3 bucket are used for running DeepVariant. 


## INPUT PARAMETERS

This pipeline does not change anyting in the original DeepVariant. It is just translated in nextflow. 
For this reason, if you are new to the pipeline and interesten in how it works and which are the prerequisites to run it, please refer to the google offical DeepVarint documentation.

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
    
    

    
    
