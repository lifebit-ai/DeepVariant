#for all bam files run deepvariant ( once at a time! No parallelization built!)
for  bam in /dv2/input/*bam
do
  start=$(date +%s.%N)
  LOGDIR=./logs
  N_SHARDS=$1
  OUTPUT_DIR=quickstart-output
  mkdir -p "${OUTPUT_DIR}"

  mkdir -p "${LOGDIR}"
  time seq 0 $((N_SHARDS-1)) | \
    parallel --eta --halt 2 --joblog "${LOGDIR}/log" --res "${LOGDIR}" \
    python ./opt/deepvariant/bin/make_examples.zip  \
      --mode calling \
      --ref /dv2/input/hg19.fa \
      --reads $bam \
      --examples "${OUTPUT_DIR}/examples.tfrecord@${N_SHARDS}.gz" \
      --task {}

      python /opt/deepvariant/bin/call_variants.zip \
       --outfile "quickstart-output/call_variants_output.tfrecord.gz" \
       --examples "quickstart-output/examples.tfrecord@${N_SHARDS}.gz" \
       --checkpoint /dv2/models/model.ckpt


      python /opt/deepvariant/bin/postprocess_variants.zip \
         --ref /dv2/input/hg19.fa.gz \
         --infile "quickstart-output/call_variants_output.tfrecord.gz" \
         --outfile "quickstart-output/output.vcf"


  end=$(date +%s.%N)
  runtime=$(python -c "print(${end} - ${start})")
  echo "$bam, Runtime was $runtime" >> $PWD/trace.txt

done

