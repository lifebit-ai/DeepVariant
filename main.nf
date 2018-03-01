/*
*  Google DeepVariant as a Nextflow pipeline!
*
*  LifeBit Biotech, 2018.
*
*/

import java.util.List;

/*
* INPUT FOLDER
* example of content:
* model.ckpt.data-00000-of-00001	model.ckpt.index		model.ckpt.meta
*/
params.modelFolder="$baseDir/data/models";
//params.modelFolder="s3://deepvariant-test/models"
params.modelName="model.ckpt";

// Names of the file to be usedcd
params.bam_definition=".bam";
// If needed
params.ref_name="ucsc.hg19.chr20.unittest.fasta";

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
model=file("${params.modelFolder}");




/*
* Generating and/or parsing the bai, fai, fai.gz, ... needed for the makeExamples process.
* Preprocessing bam files adding sample info if not contained
*/

params.bam_folder="/Users/luisasantus/DeepVariant/data/input";
params.bai_folder="none";
params.fasta="/Users/luisasantus/DeepVariant/data/input/ucsc.hg19.chr20.unittest.fasta";
params.fai="none";

fasta=file(params.fasta);

bamChannel= Channel.fromPath("${params.bam_folder}/*.bam");
bai_folder=Channel.fromPath("${params.bai_folder}/*.bai");



process preprocessFASTA{

  container 'luisas/samtools'

  input:
  file fasta from fasta

  output:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi") into fastaChannel

  script:
  """
  ## if not given --> cretae index for fasta file (fai)ls
  [[ "${params.fai}" == "none" ]] && samtools faidx $fasta && bgzip -c -i ${fasta} > ${fasta}.gz && samtools faidx "${fasta}.gz"  || cp ${params.fai} .
  """

}


process preprocessBAM{

  container 'luisas/samtools'
  publishDir "$baseDir/out"


  input:
  file bam from bamChannel

  output:
  set file(bam), file("${bam}.bai") into completeChannel

  script:
  """
    ## if not bam files
    [[ "${params.bai_folder}" == "none" ]] && samtools index $bam

    [[ `samtools view -H $bam | grep '@RG' | wc -l`   > 0 ]] || java -jar picard.jar AddOrReplaceReadGroups \
      I=$bam \
      O=$bam \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
  """

}



/*
* Getting bam files and converting them to images ( named examples )
*/
process makeExamples{

  input:
  set file(bam), file("${bam}.bai") from completeChannel
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi") from fastaChannel


  output:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"),val(bam), file("shardedExamples") into examples

  script:
  """
    mkdir shardedExamples

    time seq 0 $numberShardsMinusOne | \
    parallel --eta --halt 2 \
      python /opt/deepvariant/bin/make_examples.zip \
      --mode calling \
      --ref "${fasta} "\
      --reads "$bam" \
      --regions "${params.regions}" \
      --examples "shardedExamples/examples.tfrecord@${params.n_shards}.gz"\
      --task {}
  """
}

/*
* Doing the variant calling based on the ML trained model.
*/

process call_variants{

  input:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"),val(bam), file("shardedExamples") from examples
  file 'dv2/models' from model

  output:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"), val(bam), file('call_variants_output.tfrecord') into called_variants

  script:
  """
  /opt/deepvariant/bin/call_variants \
    --outfile call_variants_output.tfrecord \
    --examples shardedExamples/examples.tfrecord@${params.n_shards}.gz \
    --checkpoint dv2/models/${params.modelName}
  """

}

/*
* Trasforming the variant calling output (tfrecord file) into a standard vcf file.
*/
process postprocess_variants{

  publishDir params.resultdir, mode: 'copy'

  input:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"), val(bam),file('call_variants_output.tfrecord') from called_variants

  output:
   set val(bam),file("$bam--${params.ref_name}.vcf") into postout

  script:
  """
    /opt/deepvariant/bin/postprocess_variants \
    --ref "${fasta}.gz" \
    --infile call_variants_output.tfrecord \
    --outfile "$bam--${params.ref_name}.vcf"
  """
}


workflow.onComplete {
    println ( workflow.success ? "Done! \nYou can find your results in $baseDir/${params.resultdir}" : "Oops .. something went wrong" )
}
