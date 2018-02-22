/*
*  Google DeepVariant as a Nextflow pipeline!
*
*  LifeBit Biotech, 2018.
*
*/



/*
* INPUT FOLDER
*
* example of content:
*
* NA12878_S1.chr20.10_10p1mb.bam			ucsc.hg19.chr20.unittest.fasta
* NA12878_S1.chr20.10_10p1mb.bam.bai		ucsc.hg19.chr20.unittest.fasta.fai
* test_nist.b37_chr20_100kbp_at_10mb.bed		ucsc.hg19.chr20.unittest.fasta.gz
* test_nist.b37_chr20_100kbp_at_10mb.vcf.gz	ucsc.hg19.chr20.unittest.fasta.gz.fai
* test_nist.b37_chr20_100kbp_at_10mb.vcf.gz.tbi	ucsc.hg19.chr20.unittest.fasta.gz.gzi
*/

//params.input="$baseDir/data/input"
params.input="s3://deepvariant-test/input"

/*
* INPUT FOLDER
*
* example of content:
* model.ckpt.data-00000-of-00001	model.ckpt.index		model.ckpt.meta
*/
//params.modelFolder="$baseDir/data/model";
params.modelFolder="s3://deepvariant-test/models"
params.modelName="model.ckpt";

// Names of the file to be used
params.bam_name="NA12878_S1.chr20.10_10p1mb.bam";
params.ref_name="ucsc.hg19.chr20.unittest.fasta";

//OTHER PARAMETERS
params.regions="chr20:10,000,000-10,010,000";



params.n_shards="1";
params.nameOutput="output";
//Name of the directory in which the vcf result will be stored
params.resultdir = "RESULTS-DeepVariant";
params.outputdir="quickstart-output"
params.log="./logs"

//Files
input_folder= file(params.input);
model=file("${params.modelFolder}");




process makeExamples{

  input:
  file 'dv2/input' from input_folder

  output:
  file 'examples.tfrecord.gz' into examples

  script:
  """
  /opt/deepvariant/bin/make_examples \
  --mode calling   \
  --ref "dv2/input/${params.ref_name}"   \
  --reads "dv2/input/${params.bam_name}" \
  --regions "${params.regions}" \
  --examples "examples.tfrecord.gz"
  """
}


process call_variants{

  input:
  file 'examples.tfrecord.gz' from examples
  file 'dv2/models' from model
  file 'dv2/input' from input_folder

  output:
  file 'call_variants_output.tfrecord' into called_variants

  script:
  """
  /opt/deepvariant/bin/call_variants \
    --outfile call_variants_output.tfrecord \
    --examples examples.tfrecord.gz \
    --checkpoint dv2/models/${params.modelName}
  """

}

process postprocess_variants{

  publishDir params.resultdir, mode: 'copy'

  input:
    file 'dv2/input' from input_folder
    file 'call_variants_output.tfrecord' from called_variants

  output:
   file 'output.vcf' into output

  script:
  """
    /opt/deepvariant/bin/postprocess_variants \
    --ref 'dv2/input/${params.ref_name}.gz' \
    --infile call_variants_output.tfrecord \
    --outfile output.vcf
  """
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Something went wrong" )
}
