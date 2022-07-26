FROM rocker/tidyverse:4.1.2
LABEL maintainer="danie"
WORKDIR /home/root
RUN apt-get update && apt-get install -y \
  htop \
  libbz2-dev \
  libcairo2-dev \
  libcurl4-openssl-dev \
  libfreetype6-dev \
  libfribidi-dev \
  libgsl-dev \
  libharfbuzz-dev \
  libjpeg-dev \
  liblzma-dev \
  libnode-dev\
  libpixman-1-dev \
  libpng-dev \
  libproj-dev \
  librsvg2-dev \
  libtiff5-dev \
  libv8-dev\
  libx11-dev \
  libxt-dev \
  libz-dev \
  tmux \
  zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*


RUN Rscript -e "install.packages(\"BiocManager\")"
RUN Rscript -e "BiocManager::install(c( \
    \"Rhtslib\", \
    \"Rsamtools\", \
    \"GenomicAlignments\", \
    \"rtracklayer\", \
    \"BSgenome\", \
    \"GenomicFeatures\", \
    \"BiocGenerics\", \
    \"MatrixGenerics\", \
    \"GenomeInfoDb\", \
    \"GenomicRanges\", \
    \"SummarizedExperiment\", \
    \"VariantAnnotation\" \
))"


##put the important packages on the PATH

ENV biocmanager_home /home/root/BiocManager/
ENV PATH ${biocmanager_home}/bin:$PATH

ENV variantannotation_home /home/root/VariantAnnotation/
ENV PATH ${variantannotation_home}/bin:$PATH


##need to change the script not to use the table, but the main.nf file instead.
##the script will be addressed later in the main.nf file.

