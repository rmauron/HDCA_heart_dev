FROM rocker/rstudio:4.3.1

RUN apt-get update && apt-get install -y \
	libproj-dev \
	libgdal-dev \
	libmagick++-dev \
	libglpk-dev \
	libxml2 \
	libxt6 \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	libpcre3-dev \
	libicu-dev \
	libjpeg-dev \
	libpng-dev \
	libxml2-dev \
	libglpk-dev \
	&& rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Copy project files
RUN mkdir -p /home/rstudio/project
WORKDIR /home/rstudio/project
COPY renv.lock renv.lock

# Change ownership
RUN chown -R rstudio:rstudio /home/rstudio/

# Install renv and restore
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "library(renv)"

# Manually install Bioconductor packages
RUN R -e "BiocManager::install(ask = FALSE)"
RUN R -e "BiocManager::install(c('S4Arrays', 'DelayedArray', 'AUCell'), ask = FALSE)"

# Mannually install SCENIC
RUN R -e "install.packages(c('remotes'), repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/ggrepel/ggrepel_0.3.tar.gz', repos = NULL, type = 'source')"
RUN R -e 'remotes::install_github("aertslab/SCENIC")'

# Restore environment
RUN R -e "renv::restore()"