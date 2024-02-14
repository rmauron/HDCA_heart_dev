# HDCA_heart_dev

This repository contains the scripts used for the analysis of the single-cell RNA seq data of the embryonic heart developmental atlas (HDCA). 
Since several environments have been used along the project, below you will find a workflow chart presenting the pipeline used with the corresponding scripts.

## Workflow chart
On the diagram, find the workflow used for the analysis with each major step presented as the title of the box, the path were the corresponding script are located and the color representing the corresponding environment.

- All the scripts are found in the [code folder](./code).
- All the environments are found in the [environments folder](./environments).

![HDCA_heart_pipeline drawio](https://github.com/rmauron/HDCA_heart_dev/assets/92672952/ed6a7361-ee9c-4c09-8dfc-fb71b57342ed)


## Data
- The preprocessed data, the intermediate .Rds object and metadata required to reproduce the analysis are found on Mendeley Data at [this DOI](doi: 10.17632/fhtb99mdzd.1).
- The raw sequencing data are availbale at the European Genome-Phenome Archive and can be accessed upon request at [this EGA](https://ega-archive.org/).
- Metadata as annotation files or stereoscope output are also shared in the [metadata folder](./metadata).

## Resources
The analyses presented in this repository were predominantly conducted on a MacBook Pro M2 Max chip (2023), 32 GB of memory, running Ventura 13.0; however, some computations were also executed on a private server for enhanced performance and scalability.

Allthough an extensive effort was attributed to reproducibility, some system dependencies might slightly affect the results.

## Citation
TBA
