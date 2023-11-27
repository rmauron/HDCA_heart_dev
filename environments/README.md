# Environments

The environments used in this project can be restored throught the YAML files for the conda environments and with the docker image found on Dockerhub for the container.

## Conda Environments

### Restore conda environments:

1. Download the `.yml` files (e.g. after git clone the project)
2. Open a terminal where the `.yml` files are downloaded
3. Run:

    ```
    conda env create -f environment.yml
    ```

4. Change `environment.yml` with the names of the `.yml` to restore (`r-semla.yml`, `scFates.yml`, `scVelo.yml`, `stereoscope.yml`)


### Use conda environemnts:

Use the conda environments by running:

   ```bash
   conda activate environment
   ```



## Docker Container

### Building the Docker Container

Most of the analysis is run on a docker container. The corresponding image is found on (dockerhub)[https://hub.docker.com/].


1. **Pull the Docker Image from Docker Hub:**

   ```bash
   docker pull yourusername/your-repo-name:tag
   ```
