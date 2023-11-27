# Environments

The environments used in this project can be restaured throught the YAML files for the conda environments and with the docker image found on Dockerhub for the container.

## Conda Environments

### Restaure conda environments:
Restaure the conda environments by:

1. Download the `.yml` files (e.g. after git clone the project)
2. Open a terminal where the `.yml` files are downloaded
3. Run

    ```
    conda env create -f environment.yml
    ```

Change `environment.yml` with the names of the `.yml` to restaure (`r-semla.yml`, `scFates.yml`, `scVelo.yml`, `stereoscope.yml`)

### Use conda environemnts:

Use the conda environments by running:

    ```
    conda env create -f environment.yml
    ```



## Docker Container

### Building the Docker Container

You can also use Docker to run your project in a containerized environment.

1. **Pull the Docker Image from Docker Hub:**

   ```bash
   docker pull yourusername/your-repo-name:tag
   ```
