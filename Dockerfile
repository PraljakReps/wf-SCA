FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:6839-main

RUN apt-get install -y curl unzip
# Or use managed library distributions through the container OS's package
# manager.
RUN apt-get update -y && \
    apt-get install -y git

RUN apt-get install gcc

RUN pip install scipy matplotlib biopython

# install FASTA from GITHUB
RUN apt-get -y install emboss

RUN git clone https://github.com/PraljakReps/pySCA_temp.git
RUN git clone https://gitlab.com/ranganathanlab/pySCA-data.git
RUN cp -r ./pySCA-data/* ./pySCA_temp/data
RUN cp -r ./pySCA_temp/pysca ./pySCA_temp/bin/pysca



# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN python3 -m pip install --upgrade latch
WORKDIR /root
