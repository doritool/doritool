FROM willmclaren/ensembl-vep
MAINTAINER "Miguel Madrid" mmadrid@cnio.es

# https://github.com/Ensembl/ensembl-vep/blob/release/88/docker/Dockerfile
# https://hub.docker.com/r/willmclaren/ensembl-vep/
# https://github.com/rocker-org/rocker/blob/master/r-base/Dockerfile
# https://hub.docker.com/_/r-base/

USER root

# Minimun packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    ed less	locales wget ca-certificates fonts-texgyre \
    libcurl4-gnutls-dev libxml2-dev

# Config locale
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# Install R and BioConductor packages
RUN apt-get install -y --no-install-recommends \
      littler \
      r-cran-littler \
      r-base \
      r-base-dev \
      r-recommended \
      && echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
      && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
      && ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
      && ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
      && ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
      && ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
      && install.r docopt \
      && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
      && rm -rf /var/lib/apt/lists/* \
      && echo 'install.packages(c("RCurl", "XML", "gProfileR"), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R \
      && echo 'source("https://bioconductor.org/biocLite.R")' >> /tmp/packages.R \
      && echo 'biocLite(c("org.Hs.eg.db","FGNet","AnnotationDbi","topGO","KEGGprofile"))' >> /tmp/packages.R \
      && Rscript /tmp/packages.R \
      && rm -rf /tmp/*/downloaded_packages/ /tmp/packages.R \
      && rm -rf /var/lib/apt/lists/*

USER vep

RUN mkdir -p ${HOME}/.vep/Plugins/config/ && cd ${HOME}/.vep/Plugins/config && \
    cd /tmp/ && git clone https://github.com/Ensembl/VEP_plugins.git && \
    cp -r VEP_plugins/* ${HOME}/.vep/Plugins/

WORKDIR $HOME/doritool

ENTRYPOINT ["./code/DoriTool.sh"]

