FROM ubuntu:latest


#intall python2.7
RUN  apt update &&  apt install -y python2.7 python-pip

RUN apt install -y wget

#install gsutils
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-183.0.0-linux-x86_64.tar.gz &&\
tar xf google-cloud-sdk-183.0.0-linux-x86_64.tar.gz &&\
./google-cloud-sdk/install.sh

RUN apt-get install -y git
#get repo
RUN git clone https://github.com/google/deepvariant



#prepare for deepVariant
RUN cd deepvariant && bash ./build-prereq.sh
RUN bash ./build_and_test.sh


#prep


ENV BUCKET gs://deepvariant
ENV BIN_VERSION 0.5.1
ENV MODEL_VERSION 0.5.0
ENV MODEL_CL 182548131

# Note that we don't specify the CL number for the binary, only the bin version.
ENV BIN_BUCKET ${BUCKET}/binaries/DeepVariant/${BIN_VERSION}/DeepVariant-${BIN_VERSION}+cl-*
ENV MODEL_NAME DeepVariant-inception_v3-${MODEL_VERSION}+cl-${MODEL_CL}.data-wgs_standard
ENV MODEL_BUCKET ${BUCKET}/models/DeepVariant/${MODEL_VERSION}/${MODEL_NAME}
ENV DATA_BUCKET ${BUCKET}/quickstart-testdata

RUN mkdir -p bin && gsutil -m cp "${BIN_BUCKET}/*" bin/ && chmod a+x bin/*
