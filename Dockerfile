FROM satijalab/seurat

#RUN useradd -ms /bin/bash gpuser
USER root
RUN mkdir /opt/genepatt
WORKDIR /opt/genepatt

RUN mkdir /src
RUN mkdir /testdata

#USER gpuser
RUN Rscript  -e "install.packages('devtools')"
RUN Rscript  -e "install.packages('harmony')"
RUN Rscript  -e "install.packages('optparse')"
RUN Rscript  -e "install.packages('xfun')"
RUN Rscript  -e "install.packages('dplyr')"

COPY src/* /opt/genepatt/src/
