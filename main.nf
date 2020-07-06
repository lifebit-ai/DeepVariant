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
params.exome="";
if(params.exome){
  model=file("s3://deepvariant-data/models/exome");
}
else{
  model=file("${params.modelFolder}");
}


/*--------------------------------------------------
  Using the BED file
---------------------------------------------------*/
params.bed=""
if(params.exome){
  assert (params.bed != true) && (params.bed != null) : "please specify --bed option (--bed bedfile)"
  bedfile=file("${params.bed}")
}

/*--------------------------------------------------
  Cores of the machine --> used for process makeExamples
  default:2
---------------------------------------------------*/
int cores = Runtime.getRuntime().availableProcessors();
params.j=cores
numberShardsMinusOne=params.j-1;

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

params.hg19="true";
params.h38="";
params.test="";
params.hg19chr20="";
params.grch37primary="";
params.hs37d5="";

params.fasta="nofasta";
params.fai="nofai";
params.fastagz="nofastagz";
params.gzfai="nogzfai";
params.gzi="nogzi";

if(!("nofasta").equals(params.fasta)){
  fasta=file(params.fasta)
  fai=file(params.fai);
  fastagz=file(params.fastagz);
  gzfai=file(params.gzfai);
  gzi=file(params.gzi);
}
else if(params.h38 ){
  fasta=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa");
  fai=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.fai");
  fastagz=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.gz");
  gzfai=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.gz.fai");
  gzi=file("s3://deepvariant-data/genomes/h38/GRCh38.p10.genome.fa.gz.gzi");
}
else if(params.test){
  fasta=file("$baseDir/testdata/ucsc.hg19.chr20.unittest.fasta");
  fai=file("$baseDir/testdata/ucsc.hg19.chr20.unittest.fasta.fai");
  fastagz=file("$baseDir/testdata/ucsc.hg19.chr20.unittest.fasta.gz");
  gzfai=file("$baseDir/testdata/ucsc.hg19.chr20.unittest.fasta.gz.fai");
  gzi=file("$baseDir/testdata/ucsc.hg19.chr20.unittest.fasta.gz.gzi");
}
else if(params.hs37d5){
  fasta=file("s3://deepvariant-data/genomes/hs37d5/hs37d5.fa");
  fai=file("s3://deepvariant-data/genomes/hs37d5/hs37d5.fa.fai");
  fastagz=file("s3://deepvariant-data/genomes/hs37d5/hs37d5.fa.gz");
  gzfai=file("s3://deepvariant-data/genomes/hs37d5/hs37d5.fa.gz.fai");
  gzi=file("s3://deepvariant-data/genomes/hs37d5/hs37d5.fa.gz.gzi");
}
else if(params.grch37primary){
  fasta=file("s3://deepvariant-data/genomes/GRCh37.dna.primary/Homo_sapiens.GRCh37.dna.primary_assembly.fa");
  fai=file("s3://deepvariant-data/genomes/GRCh37.dna.primary/Homo_sapiens.GRCh37.dna.primary_assembly.fa.fai");
  fastagz=file("s3://deepvariant-data/genomes/GRCh37.dna.primary/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz");
  gzfai=file("s3://deepvariant-data/genomes/GRCh37.dna.primary/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz.fai");
  gzi=file("s3://deepvariant-data/genomes/GRCh37.dna.primary/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz.gzi");
}
else if(params.hg19chr20 ){
  fasta=file("s3://deepvariant-data/genomes/hg19chr20/chr20.fa");
  fai=file("s3://deepvariant-data/genomes/hg19chr20/chr20.fa.fai");
  fastagz=file("s3://deepvariant-data/genomes/hg19chr20/chr20.fa.gz");
  gzfai=file("s3://deepvariant-data/genomes/hg19chr20/chr20.fa.gz.fai");
  gzi=file("s3://deepvariant-data/genomes/hg19chr20/chr20.fa.gz.gzi");
}
else if(params.hg19 ){
  fasta=file("s3://deepvariant-data/genomes/hg19/hg19.fa");
  fai=file("s3://deepvariant-data/genomes/hg19/hg19.fa.fai");
  fastagz=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz");
  gzfai=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz.fai");
  gzi=file("s3://deepvariant-data/genomes/hg19/hg19.fa.gz.gzi");
}

else{
  System.out.println(" --fasta \"/path/to/your/genome\"  params is required and was not found! ");
  System.out.println(" or you can use standard genome versions by typing --hg19 or --h38 ");
  System.exit(0);
}



/*--------------------------------------------------
  Bam related input files
---------------------------------------------------*/

params.getBai="false";

if(params.test){
    params.bam_folder="$baseDir/testdata"
}

assert (params.bam_folder != true) && (params.bam_folder != null) : "please specify --bam_folder option (--bam_folder bamfolder)"


params.bam_file_prefix="*"

if( !("false").equals(params.getBai)){
  Channel.fromFilePairs("${params.bam_folder}/${params.bam_file_prefix}*.{bam,bam.bai}").set{bamChannel}
}else{
  Channel.fromPath("${params.bam_folder}/${params.bam_file_prefix}*.bam").map{ file -> tuple(file.name, file) }.set{bamChannel}
}

/*--------------------------------------------------
  Output directory
---------------------------------------------------*/
params.resultdir = "results";

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

  container 'lifebitai/preprocessingvctools'
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
  [[ "${params.fai}"=="nofai" ]] &&  samtools faidx $fasta || echo " fai file of user is used, not created"
  [[ "${params.fastagz}"=="nofastagz" ]]  && bgzip -c ${fasta} > ${fasta}.gz || echo "fasta.gz file of user is used, not created "
  [[ "${params.gzi}"=="nogzi" ]] && bgzip -c -i ${fasta} > ${fasta}.gz || echo "gzi file of user is used, not created"
  [[ "${params.gzfai}"=="nogzfai" ]] && samtools faidx "${fasta}.gz" || echo "gz.fai file of user is used, not created"
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
  container 'lifebitai/samtools'
  publishDir "$baseDir/sampleDerivatives"

  input:
  set val(prefix), file(bam) from bamChannel
  output:
  set file("ready/${bam[0]}"), file("ready/${bam[0]}.bai") into completeChannel, completeStats
  script:
  """
	  mkdir ready
  [[ `samtools view -H ${bam[0]} | grep '@RG' | wc -l`   > 0 ]] && { mv $bam ready;}|| { picard AddOrReplaceReadGroups \
    I=${bam[0]} \
    O=ready/${bam[0]} \
    RGID=${params.rgid} \
    RGLB=${params.rglb} \
    RGPL=${params.rgpl} \
    RGPU=${params.rgpu} \
    RGSM=${params.rgsm};}
    cd ready ;samtools index ${bam[0]};
  """
}


process BAMstats{

  tag "${bam[0]}"
  container 'lifebitai/samtools'

  input:
  set file(bam), file(bai) from completeStats
  output:
  file("*") into bam_multiqc
  script:
  """
  samtools stats $bam > stats.txt
  samtools flagstat $bam > flagstat.txt
  samtools idxstats $bam > idxstats.txt
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

if(params.bed){
  process makeExamples_with_bed{

      tag "${bam[1]}"
    cpus params.j

    input:
      set file(fasta), file(bam) from all_fa_bam
      file bedfile from bedfile
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
        --examples shardedExamples/examples.tfrecord@!{params.j}.gz\
        --regions !{bedfile} \
        --task {}
    '''
  }
}
else{
  process makeExamples{

    tag "${bam[1]}"
    cpus params.j

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
        --examples shardedExamples/examples.tfrecord@!{params.j}.gz\
        --task {}
    '''
  }
}
/********************************************************************
  process call_variants
  Doing the variant calling based on the ML trained model.
********************************************************************/



process call_variants{


  tag "${bam}"
  cpus params.j

  input:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"),val(bam), file("shardedExamples") from examples
  file 'dv2/models' from model
  output:
  set file(fasta),file("${fasta}.fai"),file("${fasta}.gz"),file("${fasta}.gz.fai"), file("${fasta}.gz.gzi"), val(bam), file('call_variants_output.tfrecord') into called_variants
  script:
  """
  /opt/deepvariant/bin/call_variants \
    --outfile call_variants_output.tfrecord \
    --examples shardedExamples/examples.tfrecord@${params.j}.gz \
    --checkpoint dv2/models/${params.modelName} \
    --num_readers ${params.j}
  """
}



/********************************************************************
  process call_variants
  Trasforming the variant calling output (tfrecord file) into a standard vcf file.
********************************************************************/

process postprocess_variants{


  tag "$bam"
  cpus params.j

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


process vcftools{
  tag "$vcf"

  container 'lifebitai/vcftools:latest'

  input:
  set val(bam),file(vcf) from postout
  output:
  file("*") into vcfout

  script:
  """
  vcftools --vcf $vcf --TsTv-summary
  vcftools --vcf $vcf --TsTv-by-count
  vcftools --vcf $vcf --TsTv-by-qual
  # remove rows containing 'inf' which breaks multiqc report
  sed -i '/inf/d' out.TsTv.qual
  """
}

process multiqc{
  tag "multiqc_report.html"

  publishDir "${params.resultdir}/MultiQC", mode: 'copy'
  container 'lifebitai/multiqc:v1.7'

  input:
  file(vcfout) from vcfout
  file(bamout) from bam_multiqc
  output:
  file("*") into multiqc

  script:
  """
  multiqc . -m vcftools -m samtools
  """
}


workflow.onComplete {
    println ( workflow.success ? "Done! \nYou can find your results in $baseDir/${params.resultdir}" : "Oops .. something went wrong" )
}
