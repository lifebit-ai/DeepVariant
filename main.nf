/*
*  Google DeepVariant as a Nextflow pipeline!
*
*  LifeBit Biotech, 2018.
*
*/

//Name of the directory in which the vcf result will be stored
params.resultdir = "RESULTS-DeepVariant";

params.regions="chr20:10,000,000-10,010,000";
params.n_shards="4";

params.input="/Users/luisasantus/awsProva/input"
input_file= file(params.input);

params.bam_name="NA12878_S1.chr20.10_10p1mb.bam";
params.ref_name="ucsc.hg19.chr20.unittest.fasta";


params.bam="/Users/luisasantus/awsProva/data/NA12878_S1.chr20.10_10p1mb.bam";
params.bai="/Users/luisasantus/awsProva/data/NA12878_S1.chr20.10_10p1mb.bam.bai";
params.fai="/Users/luisasantus/awsProva/data/ucsc.hg19.chr20.unittest.fasta.fai"
params.ref="/Users/luisasantus/awsProva/data/ucsc.hg19.chr20.unittest.fasta";
params.ref_gz="/Users/luisasantus/awsProva/data/ucsc.hg19.chr20.unittest.fasta.gz";
params.ref_gz_gzi="/Users/luisasantus/awsProva/data/ucsc.hg19.chr20.unittest.fasta.gz.gzi";

params.ref_gz_fai="/Users/luisasantus/awsProva/data/ucsc.hg19.chr20.unittest.fasta.gz.fai";

params.modelFolder="/Users/luisasantus/awsProva/models";
model=file("${params.modelFolder}");

model1=("${params.modelFolder}/model.ckpt.data-00000-of-00001");
model2=("${params.modelFolder}/model.ckpt.index");
model3=("${params.modelFolder}/model.ckpt.meta");

params.modelName="model.ckpt";

params.output_dir ="quickstart-output"





ref_gz_file = file(params.ref_gz);
ref_gz_gzi_file = file(params.ref_gz_gzi);
ref_gz_fai_file = file(params.ref_gz_fai);

ref_channel = Channel.from(file(params.ref))
bam_channel = Channel.from(file(params.bam))
bai_channel = Channel.from(file(params.bai))
fai_channel = Channel.from(file(params.fai))


process makeExamples{

  input:
  file 'dv2/input' from input_file
  //file 'file.bam' from bam_channel
  //file 'file.bai' from bai_channel
  //file 'reference.fasta.fai' from fai_channel
  //file 'reference.fasta' from ref_channel

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
  file 'dv2/input' from input_file

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
    file 'reference.fasta.gz' from ref_gz_file
    file 'reference.fasta.gz.gzi' from ref_gz_gzi_file
    file 'reference.fasta.gz.fai' from ref_gz_fai_file
    file 'call_variants_output.tfrecord' from called_variants

  output:
   file 'output.vcf' into output

  script:
  """
    /opt/deepvariant/bin/postprocess_variants \
    --ref 'reference.fasta.gz' \
    --infile call_variants_output.tfrecord \
    --outfile output.vcf
  """
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
