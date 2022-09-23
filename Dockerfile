FROM rocker/r-ubuntu:20.04

RUN useradd -ms /bin/bash gpuser
USER root
RUN mkdir /opt/genepatt
WORKDIR /opt/genepatt

USER root
RUN mkdir /src
RUN mkdir /testdata

#USER gpuser
RUN Rscript  -e "install.packages('devtools')"
RUN Rscript  -e "install.packages('harmony')"
RUN Rscript  -e "install.packages('optparse')"
RUN Rscript  -e "install.packages('xfun')"
RUN Rscript  -e "install.packages('png')"
RUN Rscript  -e "install.packages('Seurat')"
RUN Rscript  -e "install.packages('dplyr')"

COPY src/* /opt/genepatt/src/
