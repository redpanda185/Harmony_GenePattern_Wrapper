FROM satijalab/seurat:4.1.0

#RUN useradd -ms /bin/bash gpuser
USER root
RUN mkdir /opt/genepatt
WORKDIR /opt/genepatt

RUN mkdir /src
RUN mkdir /testdata

#USER gpuser
RUN Rscript  -e "install.packages('devtools', version = '2.4.5')"
RUN Rscript  -e "install.packages('harmony', version = '0.1.1')"
RUN Rscript  -e "install.packages('optparse', version = '1.7.3')"
RUN Rscript  -e "install.packages('xfun', version = '0.36')"
RUN Rscript  -e "install.packages('dplyr', version = '1.0.10')"

COPY src/* /opt/genepatt/src/
