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
* ucsc.hg19.chr20.unittest.fasta.gz
* ucsc.hg19.chr20.unittest.fasta.gz.fai
* ucsc.hg19.chr20.unittest.fasta.gz.gzi
*/

params.input="$baseDir/data/input"
//params.input="s3://deepvariant-test/input"

/*
* INPUT FOLDER
*
* example of content:
* model.ckpt.data-00000-of-00001	model.ckpt.index		model.ckpt.meta
*/
params.modelFolder="$baseDir/data/models";
//params.modelFolder="s3://deepvariant-test/models"
params.modelName="model.ckpt";

// Names of the file to be used
//params.bam_definition="NA12878_S1.chr20.10_10p1mb.bam";
params.bam_definition=".bam";
// If needed
//params.bam_definition=".bam";
params.ref_name="ucsc.hg19.chr20.unittest.fasta";

//Obtain all bam files in input file directory
bamchannel = Channel.create()
File dir = new File("${params.input}");
for(File f : dir.listFiles()){
  if(f.getName().endsWith(params.bam_definition)){
        bamchannel<<f.getName();
  }
}
bamchannel.close();

//OTHER PARAMETERS
params.regions="chr20:10,000,000-10,010,000";

//Number of cores in the Machine
params.n_shards=1;
numberShardsMinusOne=params.n_shards-1;
shardsChannel= Channel.from( 0..params.n_shards);


//Name of the directory in which the vcf result will be stored
params.nameOutput="${params.ref_name}";
params.resultdir = "RESULTS-DeepVariant";
params.outputdir="quickstart-output"
params.log="./logs"

//Files
input_folder= file(params.input);
model=file("${params.modelFolder}");

///opt/deepvariant/bin/make_examples \

process makeExamples{

  input:
  val bam from bamchannel
  file 'dv2/input' from input_folder


  output:
  set val(bam),  file("shardedExamples") into examples

  script:
  """
    mkdir shardedExamples

    time seq 0 $numberShardsMinusOne | \
    parallel --eta --halt 2 \
      python /opt/deepvariant/bin/make_examples.zip \
      --mode calling \
      --ref "dv2/input/${params.ref_name} "\
      --reads "dv2/input/$bam" \
      --regions "${params.regions}" \
      --examples "shardedExamples/examples.tfrecord@${params.n_shards}.gz"\
      --task {}
  """
}


process call_variants{

  input:
   set val(bam), file('shardedExamples') from examples
  file 'dv2/models' from model
  file 'dv2/input' from input_folder

  output:
   set val(bam), file('call_variants_output.tfrecord') into called_variants

  script:
  """
  /opt/deepvariant/bin/call_variants \
    --outfile call_variants_output.tfrecord \
    --examples shardedExamples/examples.tfrecord@${params.n_shards}.gz \
    --checkpoint dv2/models/${params.modelName}
  """

}

process postprocess_variants{

  publishDir params.resultdir, mode: 'copy'

  input:
   file 'dv2/input' from input_folder
    set val(bam),file('call_variants_output.tfrecord') from called_variants

  output:
   set val(bam),file("$bam--${params.ref_name}.vcf") into postout

  script:
  """
    /opt/deepvariant/bin/postprocess_variants \
    --ref 'dv2/input/${params.ref_name}.gz' \
    --infile call_variants_output.tfrecord \
    --outfile "$bam--${params.ref_name}.vcf"
  """
}


workflow.onComplete {
    println ( workflow.success ? "Done! \n You can find your results in $baseDir/${params.resultdir}" : "Oops .. something went wrong" )
}
