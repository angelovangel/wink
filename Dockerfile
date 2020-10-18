# all requirements for the WINK pipeline

FROM continuumio/miniconda3:4.8.2

LABEL description="Docker image containing all requirements for running WINK (nextflow part)"
LABEL maintainer="Angel Angelov <aangeloo@gmail.com>"

# the base image of continuumio does not have procps, needed by nextflow
RUN apt-get update --fix-missing && \
    apt-get install -y libz-dev procps ksh && apt-get clean -y

# make bracken binaries available
ENV PATH /Bracken:$PATH

# final stage
COPY environment.yml ./

RUN conda env update -n root -f environment.yml && conda clean -afy
