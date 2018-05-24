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
  Can be substitued with own model folder.
---------------------------------------------------*/
params.modelFolder="s3://deepvariant-data/models"
params.modelName="model.ckpt";
model=file("${params.modelFolder}");


/*--------------------------------------------------
  Cores of the machine --> used for process makeExamples
  default:2
---------------------------------------------------*/
def proc = 'lscpu -e=cpu'.execute() | 'wc -l'.execute()
proc.waitFor()
def number = "${proc.in.text}" as int
number -= 1
params.n_shards=number
numberShardsMinusOne=params.n_shards-1;
shardsChannel= Channel.from( 0..params.n_shards);


/*--------------------------------------------------
  Fasta related input files

  You can use the flag --hg19 for using the hg19 version of the Genome.
  You can use the flag --h38 for using the GRCh38.p10 version of the Genome.

  They can be passed manually, through the parameter:
  	params.fasta="/my/path/to/file";
  And if already at user's disposal:
	params.fai="/my/path/to/file";
	params.fastagz="/my/path/to/file";
	params.gzfai="/my/path/to/file";
	params.gzi="/my/path/to/file";

---------------------------------------------------*/



params.fasta="nofasta";
params.fai="nofai";
params.fastagz="nofastagz";
params.gzfai="nogzfai";
params.gzi="nogzi";


 params.hg19="";
 params.h38="";

   fasta=file("s3://deepvariant-data/genomes/hg19/hg19.fa");
   fai=file("s3://deepvariant-data/genomes/hg19/hg19.fa.fai");
   fastagz=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz");
   gzfai=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz.fai");
   gzi=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz.gzi");

// if(params.hg19){
 //  fasta=file("s3://deepvariant-data/genomes/hg19/hg19.fa");
  // fai=file("s3://deepvariant-data/genomes/hg19/hg19.fa.fai");
 //  fastagz=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz");
 //  gzfai=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz.fai");
 //  gzi=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz.gzi");
// }
// else if(params.h38){
 //  fasta=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa");
  // fai=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.fai");
  // fastagz=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.gz");
 //  gzfai=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.gz.fai");
 //  gzi=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.gz.gzi");
// }
// else{
 //  fasta=file(params.fasta)
 //  fai=file(params.fai);
  // fastagz=file(params.fastagz);
  // gzfai=file(params.gzfai);
  // gzi=file(params.gzi);
// }

// if(("nofasta").equals(params.fasta) && !params.hg19 && !params.h38 ){
 //  System.out.println(" --fasta \"/path/to/your/genome\"  params is required and was not found! ");
 //  System.exit(0);
// }


/*--------------------------------------------------
  Bam related input files
---------------------------------------------------*/

params.bam_folder="s3://deepvariant-data/test-bam/hg19-mediumMultiple/";
params.getBai="true";

if( !("false").equals(params.getBai)){
  Channel.fromFilePairs("${params.bam_folder}/*.{bam,bam.bai}").set{bamChannel}
}else{
  Channel.fromPath("${params.bam_folder}/*.bam").map{ file -> tuple(file.name, file) }.set{bamChannel}
}

/*--------------------------------------------------
  Output directory
---------------------------------------------------*/
params.resultdir = "RESULTS-DeepVariant";

/*--------------------------------------------------
  Params for the Read Group Line to be added just in
  case its needed.
  If not given, default values are used.
---------------------------------------------------*/
params.rgid=4;
params.rglb="lib1";
params.rgpl="illumina";
params.rgpu="unit1";
params.rgsm=20;



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
  [[ $fai == "nofai" ]] &&  samtools faidx $fasta || echo " fai file of user is used, not created"
  [[ $fastagz == "nofastagz" ]]  && bgzip -c ${fasta} > ${fasta}.gz || echo "fasta.gz file of user is used, not created "
  [[ $gzfai == "nogzi" ]] && bgzip -c -i ${fasta} > ${fasta}.gz || echo "gzi file of user is used, not created"
  [[ $gzi == "nogzfai" ]] && samtools faidx "${fasta}.gz" || echo "gz.fai file of user is used, not created"
  """
}


/********************************************************************
  process preprocessBAM
  If the user gives the index files for the bam files as an input, they are used
  If not they are produced in this process given only the fasta file.
  Moreover this takes care of the read group line too.
********************************************************************/


process preprocessBAM{


  tag "${bam[0]}"
  container 'luisas/samtools'
  publishDir "$baseDir/sampleDerivatives"

  input:
  set val(prefix), file(bam) from bamChannel
  output:
  set file("ready/${bam[0]}"), file("ready/${bam[0]}.bai") into completeChannel
  script:
  """
	  mkdir ready
  [[ `samtools view -H ${bam[0]} | grep '@RG' | wc -l`   > 0 ]] && { mv $bam ready; cd ready;  }|| { java -jar /picard.jar AddOrReplaceReadGroups \
    I=${bam[0]} \
    O=ready/${bam[0]} \
    RGID=${params.rgid} \
    RGLB=${params.rglb} \
    RGPL=${params.rgpl} \
    RGPU=${params.rgpu} \
    RGSM=${params.rgsm}; cd ready ;samtools index ${bam[0]}; }
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

	Can be parallelized through the params.n_shards
	( if params.n_shards >= 1 parallelization happens automatically)
      ********************************************************************/

process makeExamples{

    tag "${bam[1]}"
  cpus params.n_shards

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
      --ref !{fasta[1]}.gz\
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


  tag "${bam}"
  cpus params.n_shards
  
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
    --checkpoint dv2/models/${params.modelName} \
    --num_readers ${params.n_shards}
  """
}



/********************************************************************
  process call_variants
  Trasforming the variant calling output (tfrecord file) into a standard vcf file.
********************************************************************/

process postprocess_variants{


  tag "$bam"
  cpus 1

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

