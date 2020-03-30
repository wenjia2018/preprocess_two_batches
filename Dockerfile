FROM rocker/verse:3.6.1

# Install system packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends git zlib1g-dev libxml2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY R/ /home/rstudio/legacy/R

RUN Rscript /home/rstudio/legacy/R/init_install_packages.R

EXPOSE 8787
