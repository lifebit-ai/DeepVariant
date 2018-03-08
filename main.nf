/*
*  Google DeepVariant as a Nextflow pipeline!
*
*  LifeBit Biotech, 2018.
*
*/

import java.util.List;

/*--------------------------------------------------
  Model folder
  Content: trained model.
  For exact information refer to documentation.
---------------------------------------------------*/
params.modelFolder="s3://deepvariant-test/models"
params.modelName="model.ckpt";
model=file("${params.modelFolder}");

/*--------------------------------------------------
  Regions
---------------------------------------------------*/
params.regions="chr20:10,000,000-10,010,000";

/*--------------------------------------------------
  Cores of the machine --> used for process makeExamples
---------------------------------------------------*/
params.n_shards=2;
numberShardsMinusOne=params.n_shards-1;
shardsChannel= Channel.from( 0..params.n_shards);


/*--------------------------------------------------
  Fasta related input files
---------------------------------------------------*/

params.fasta="s3://deepvariant-test/input/ucsc.hg19.chr20.unittest.fasta";
params.fai="nofai";
params.fastagz="nofastagz";
params.gzfai="nogzfai";
params.gzi="nogzi";

fasta=file(params.fasta)
fai=file(params.fai);
fastagz=file(params.fastagz);
gzfai=file(params.gzfai);
gzi=file(params.gzi);

/*--------------------------------------------------
  Bam related input files
---------------------------------------------------*/
params.bam_folder="s3://deepvariant-test/input/";
params.getBai="false";

if( !("false").equals(params.getBai)){
  Channel.fromFilePairs("${params.bam_folder}/*.{bam,bam.bai}").set{bamChannel}
}else{
  Channel.fromPath("${params.bam_folder}/*.bam").map{ file -> tuple(file.name, file) }.set{bamChannel}
}

/*--------------------------------------------------
  Output directory
---------------------------------------------------*/
params.resultdir = "RESULTS-DeepVariant";

/********************************************************************

  process preprocessFASTA

  Collects all the files related to the reference genome, like
  .fai,.gz ...

  If the user gives them as an input, they are used
  If not they are produced in this process given only the fasta file.
********************************************************************/


process preprocessFASTA{

  container 'luisas/samtools'
  publishDir "$baseDir/sampleDerivatives"

  input:
  file fasta from fasta
  file fai from fai
  file fastagz from fastagz
  file gzfai from gzfai
  file gzi from gzi

  output:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi") into fastaChannel

  script:
  """
  [[ "${params.fai}" == "nofai" ]] &&  samtools faidx $fasta || echo " fai file of user is used, not created"
  [[ "${params.fastagz}" == "nofastagz" ]]  && bgzip -c ${fasta} > ${fasta}.gz || echo "fasta.gz file of user is used, not created "
  [[ "${params.gzi}" == "nogzi" ]] && bgzip -c -i ${fasta} > ${fasta}.gz || echo "gzi file of user is used, not created"
  [[ "${params.gzfai}" == "nogzfai" ]] && samtools faidx "${fasta}.gz" || echo "gz.fai file of user is used, not created"
  """
}


/********************************************************************

  process preprocessBAM

  If the user gives the index files for the bam files as an input, they are used
  If not they are produced in this process given only the fasta file.

  Moreover this takes care of the read group line too.
********************************************************************/


process preprocessBAM{

  container 'luisas/samtools'
  publishDir "$baseDir/sampleDerivatives"

  input:
  set val(prefix), file(bam) from bamChannel

  output:
  set file("ready/${bam[0]}"), file("ready/${bam[0]}.bai") into completeChannel


  script:
  """
    mkdir ready

    if [ `samtools view -H ${bam[0]} | grep '@RG' | wc -l`   > 0 ]; then
      mv $bam ready
      if [ "${params.getBai}" == "false" ]; then
         samtools index "ready/${bam[0]}"
      fi
    else
      java -jar /picard.jar AddOrReplaceReadGroups \
        I="${bam[0]}" \
        O="ready/${bam[0]} \
        RGID=4 \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=20 ;
      samtools index "ready/${bam[0]}"

    fi
  

  """
}



fastaChannel.map{file -> tuple (1,file[0],file[1],file[2],file[3],file[4])}
            .set{all_fa};

completeChannel.map { file -> tuple(1,file[0],file[1]) }
               .set{all_bam};

all_fa.cross(all_bam)
      .set{all_fa_bam};



      /********************************************************************

        process makeExamples

        Getting bam files and converting them to images ( named examples )

      ********************************************************************/

process makeExamples{

  input:
    set file(fasta), file(bam) from all_fa_bam

  output:
    set file("${fasta[1]}"),file("${fasta[1]}.fai"),file("${fasta[1]}.gz"),file("${fasta[1]}.gz.fai"), file("${fasta[1]}.gz.gzi"),val("${bam[1]}"), file("shardedExamples") into examples

  shell:
  '''
    mkdir shardedExamples

    time seq 0 !{numberShardsMinusOne} | \
    parallel --eta --halt 2 \
      python /opt/deepvariant/bin/make_examples.zip \
      --mode calling \
      --ref !{fasta[1]}\
      --regions !{params.regions} \
      --reads !{bam[1]} \
      --examples shardedExamples/examples.tfrecord@!{params.n_shards}.gz\
      --task {}
  '''
}

/********************************************************************

  process call_variants

  Doing the variant calling based on the ML trained model.

********************************************************************/



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



/********************************************************************

  process call_variants

  Trasforming the variant calling output (tfrecord file) into a standard vcf file.

********************************************************************/

process postprocess_variants{

  publishDir params.resultdir, mode: 'copy'

  input:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"), val(bam),file('call_variants_output.tfrecord') from called_variants

  output:
   set val(bam),file("${bam}.vcf") into postout

  script:
  """
    /opt/deepvariant/bin/postprocess_variants \
    --ref "${fasta}.gz" \
    --infile call_variants_output.tfrecord \
    --outfile "${bam}.vcf"
  """
}


workflow.onComplete {
    println ( workflow.success ? "Done! \nYou can find your results in $baseDir/${params.resultdir}" : "Oops .. something went wrong" )
}
