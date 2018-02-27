/*
*  Google DeepVariant as a Nextflow pipeline!
*
*  LifeBit Biotech, 2018.
*
*/

import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.model.AmazonS3Exception;
import com.amazonaws.services.s3.model.Bucket;
import com.amazonaws.AmazonClientException;
import com.amazonaws.AmazonServiceException;
import com.amazonaws.auth.profile.ProfileCredentialsProvider;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3Client;
import com.amazonaws.services.s3.model.ListObjectsRequest;
import com.amazonaws.services.s3.model.ListObjectsV2Request;
import com.amazonaws.services.s3.model.ListObjectsV2Result;
import com.amazonaws.services.s3.model.ObjectListing;
import com.amazonaws.services.s3.model.S3ObjectSummary;
import java.util.List;

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

//params.input="$baseDir/data/input"
params.input="s3://deepvariant-test/input"

/*
* INPUT FOLDER
*
* example of content:
* model.ckpt.data-00000-of-00001	model.ckpt.index		model.ckpt.meta
*/
//params.modelFolder="$baseDir/data/models";
params.modelFolder="s3://deepvariant-test/models"
params.modelName="model.ckpt";

// Names of the file to be used
params.bam_definition=".bam";
// If needed
params.ref_name="ucsc.hg19.chr20.unittest.fasta";


// CONTROL CONTENT OF FOLDER
if(params.input.startsWith("s3")){
  bamchannel = Channel.create()
  String[] s3Info = params.input.replaceAll("s3://","").split("/");
  String bucketName = s3Info[0];
  String folderKey = s3Info[1];
  ListObjectsRequest listObjectsRequest =
                                  new ListObjectsRequest()
                                        .withBucketName(bucketName)
                                        .withPrefix(folderKey + "/");
      List<String> keys = new ArrayList<>();
      AmazonS3 s3Client = new AmazonS3Client(new ProfileCredentialsProvider());
      ObjectListing objects = s3Client.listObjects(listObjectsRequest);
      for (;;) {
          List<S3ObjectSummary> summaries = objects.getObjectSummaries();
          if (summaries.size() < 1) {
              break;
          }
          for(S3ObjectSummary s : summaries){
            keys.add(s.getKey());
          }
          objects = s3Client.listNextBatchOfObjects(objects);
      }
      for(String key : keys){
          if(key.endsWith(params.bam_definition)){
            bamchannel << key.replaceAll(folderKey,"").replaceAll("/", "");
          }
      }
      bamchannel.close();
}
else{
  //Obtain all bam files in input file directory
  if(params.aws == "false" ){
    bamchannel = Channel.create()
    File dir = new File("${params.input}");
    for(File f : dir.listFiles()){
      if(f.getName().endsWith(params.bam_definition)){
            bamchannel<<f.getName();
      }
    }
    bamchannel.close();
  }

}






//bamchanne=Channel.fromPath("${params.input}/*.bam").println();

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
    println ( workflow.success ? "Done! \nYou can find your results in $baseDir/${params.resultdir}" : "Oops .. something went wrong" )
}
