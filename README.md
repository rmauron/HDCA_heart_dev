# Spatiotemporal gene expression and cellular dynamics of the developing human heart

This repository contains the scripts used for the analysis of the single-cell RNA seq and the Spatial Transcriptomics data of the embryonic heart developmental atlas (HDCA). 
Since several environments have been used along the project, below you will find a workflow chart presenting the pipeline used with the corresponding scripts. 

<br>

<div align="right">
  <a href="https://hdca-sweden.scilifelab.se/" target="_blank">
  <img src="https://github.com/rmauron/HDCA_heart_dev/assets/92672952/5baae706-3452-49ad-a616-5cf34d768ad5" alt="HDCA" width="250">
  </a>
</div>



## Workflow chart
On the diagram, find the workflow used for the analysis with each major step presented as the title of the box, the path were the corresponding script are located and the color representing the corresponding environment.

- All the scripts are found in the [code folder](./code).
- All the environments are found in the [environments folder](./environments).

![HDCA_heart_pipeline](https://github.com/rmauron/HDCA_heart_dev/assets/92672952/5d2beacf-8f18-4474-bfc0-199f5d0c0041)


## Data
- The preprocessed data, some intermediate .Rds objects and metadata required to reproduce the analysis are found on Mendeley Data at these DOI: [Mendeley repo 1](https://data.mendeley.com/preview/fhtb99mdzd?a=27a510e3-60f7-40b9-968d-ecf1ca6b5ad1) and [Mendeley repo 2](https://data.mendeley.com/preview/w65jtfsvpr?a=2c7eb695-0a84-4bd7-98e8-e4be4e4ed831).
- The raw sequencing data are availbale at the European Genome-Phenome Archive (EGA) upon request at [EGAS50000001122](https://ega-archive.org/studies/EGAS50000001122) and [missing link](https://ega-archive.org/studies/EGAS50000001122).

## Resources
The analyses presented in this repository were predominantly conducted on a MacBook Pro M2 Max chip (2023), 32 GB of memory, running Ventura 13.0; however, some computations were also executed on a private server for enhanced performance and scalability.

Although an extensive effort was attributed to reproducibility, some system dependencies might slightly affect the results. Support on it is out of the scope of that study.
Find how to set up the docker container or the different environments in the [environments folder](./environments).

## Useful links
- [Nature Genetics](missing-link)
- [Interactive viewer](https://hdcaheart.serve.scilifelab.se/web/index.html)
- [biorXive](https://www.biorxiv.org/content/10.1101/2024.03.12.584577v3) (preprint)
- [Zenodo](https://zenodo.org/records/15912657)
- [Mendeley repo 1](https://data.mendeley.com/preview/fhtb99mdzd?a=27a510e3-60f7-40b9-968d-ecf1ca6b5ad1) (Cellranger, Spaceranger, metadata)
- [Mendeley repo 2](https://data.mendeley.com/preview/w65jtfsvpr?a=2c7eb695-0a84-4bd7-98e8-e4be4e4ed831) (R-objects)
- [EGA Visium](https://ega-archive.org/studies/EGAS50000001122) (raw sequencing data, available upon formal request)
- [EGA Single-cell](missing-link) (raw sequencing data, available upon formal request)

## Citation
TBA