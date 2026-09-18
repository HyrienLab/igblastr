ARG TAG=devel
FROM bioconductor/bioconductor_docker:${TAG}

RUN R -q -e 'install.packages("remotes")'
ARG TAG
RUN R -q -e "BiocManager::install('HyrienLab/igblastr@${TAG}')"
RUN R -q -e 'igblastr::install_igblast()'
