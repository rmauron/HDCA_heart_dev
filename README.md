# Spatial Dynamics of the Developing Human Heart

This repository contains the scripts used for the analysis of the single-cell RNA seq data of the embryonic heart developmental atlas (HDCA). 
Since several environments have been used along the project, below you will find a workflow chart presenting the pipeline used with the corresponding scripts. 

<div align= "right" style="padding-right: 30px;">
  <img src="https://github.com/rmauron/HDCA_heart_dev/assets/92672952/5baae706-3452-49ad-a616-5cf34d768ad5" alt="HDCA" width="300">
</div>



## Workflow chart
On the diagram, find the workflow used for the analysis with each major step presented as the title of the box, the path were the corresponding script are located and the color representing the corresponding environment.

- All the scripts are found in the [code folder](./code).
- All the environments are found in the [environments folder](./environments).

![HDCA_heart_pipeline drawio](https://github.com/rmauron/HDCA_heart_dev/assets/92672952/ed6a7361-ee9c-4c09-8dfc-fb71b57342ed)


## Data
- The preprocessed data, some intermediate .Rds objects and metadata required to reproduce the analysis are found on Mendeley Data at [these DOI (part 1)]([doi: 10.17632/fhtb99mdzd.1](https://data.mendeley.com/preview/fhtb99mdzd?a=27a510e3-60f7-40b9-968d-ecf1ca6b5ad1)) and [(part 2)]([doi: 10.17632/fhtb99mdzd.1](https://data.mendeley.com/preview/w65jtfsvpr?a=2c7eb695-0a84-4bd7-98e8-e4be4e4ed831)).
- The raw sequencing data will be availbale at the Federated European Genome-Phenome Archive (FEGA) upon request.

## Resources
The analyses presented in this repository were predominantly conducted on a MacBook Pro M2 Max chip (2023), 32 GB of memory, running Ventura 13.0; however, some computations were also executed on a private server for enhanced performance and scalability.

Although an extensive effort was attributed to reproducibility, some system dependencies might slightly affect the results. Support on it is out of the scope of that study.

## Citation
TBA
