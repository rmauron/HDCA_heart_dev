# Environments

The environments used in this project can be restored throught the YAML files for the conda environments and with the docker image found on Dockerhub for the container.


## Docker Container

### Building the Docker Container

Most of the analysis is ran on a docker container. The corresponding image is found on [Docker Hub](https://hub.docker.com/).
The most useful commands from the [documentation](https://docs.docker.com/language/java/run-containers/) to restore the container are the following:

1. **Pull the Docker Image from Docker Hub:**

   ```bash
   docker pull yourusername/your-repo-name:tag
   ```

2. **Build the Container:**

   ```bash
   sudo docker run -d -p 1337:8787 -p 3030:3030 --name hdca_devheart -e PASSWORD=<YOURPASSWORD> --memory=30g --mount type=bind,source="<SOURCEPATH>",target=/home/rstudio -e ROOT=TRUE hdcadevheart
   ```

**Where:**
- ```<YOURPASSWORD>``` is the password you want to use when connecting to the RStudio server.
- ```--memory``` is the RAM you want to allocate to the container.
- ```<SOURCEPATH>``` is the directory on your local computer that will be accessible form the server. Make sure the data and code are accessible from that directory.



## Conda Environments

Some analysis, especially the one ran with *semla* and the python scripts (*stereoscope*, *scVelo*, *scFates*) were run on separate conda environments.

**NB:** *stereoscope* was run on a private GPU server.

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

**NB:** Make sure your Conda is relatively updated to avoid trouble when restoring from the `.yml` files.

