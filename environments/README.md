# Environments

The environments used in this project can be restored throught the YAML files for the conda environments and with the docker image found on Dockerhub for the container.


[Docker Container](## Docker Container)

<br>

---

## Docker Container

### Building the Docker Container

Most of the analysis is ran on a docker container. The corresponding image is found on [Docker Hub](https://hub.docker.com/r/raphaelmauron/hdcadevheart) and was build from the <rocker/rstudio:4.3.1> basis.
The most useful commands from the [documentation](https://docs.docker.com/language/java/run-containers/) to restore the container are the following.

**NB:** the container is going to be mounted according to the ```<SOURCEPATH>``` provided in the later command. We recommand to clone the entire GitHub repository, and start the restauration from there:

   ```
   git clone https://github.com/rmauron/HDCA_heart_dev.git
   ```

and go in the project directory:

   ```
   cd HDCA_heart_dev
   ```

and proceed.

1. **Pull the Docker Image from Docker Hub:**

Control that Docker is installed:

   ```
   docker ps
   ```

and pull the image:

   ```bash
   docker pull raphaelmauron/hdcadevheart
   ```

This can take several minutes since the image containes ```R version 4.3.1``` and the required packages.

2. **Build the Container:**

   ```bash
   sudo docker run -d -p 1337:8787 -p 3030:3030 --name hdcadevheart -e PASSWORD=<YOURPASSWORD> --memory=30g --mount type=bind,source="$(pwd)",target=/home/rstudio -e ROOT=TRUE raphaelmauron/hdcadevheart:latest
   ```

**Where:**
- ```<YOURPASSWORD>``` is the password you want to use when connecting to the RStudio server.
- ```--memory``` is the RAM you want to allocate to the container.
- ```<SOURCEPATH>``` is the directory on your local computer that will be accessible form the server. Make sure the data and code are accessible from that directory. If you hare located in ```~/HDCA_heart_dev/``` you can use ```"$(pwd)"``` as ```<SOURCEPATH>```

3. **Run the Container on your web browser**

Check that the container is running:

   ```
   docker container ls -a
   ```

Access the container on your web browser with ```localhost:1337```.

Enter the username and password:

- ```Username```: rstudio
- ```Password```: what you passed as ```<YOURPASSWORD>``` on point 2 

<br>

---

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

