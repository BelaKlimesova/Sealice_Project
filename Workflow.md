Workflow
========

The steps to set up and use RStudio within a Docker container. By using Docker, we ensure a consistent and reproducible development environment, which simplifies package management and dependency resolution. This workflow will guide you through pulling the RStudio Docker image, running it, installing necessary R packages, and accessing the RStudio interface.

## 1. Run RStudio in a Docker container

```bash
docker pull rocker/rstudio
docker run -d -p 8787:8787 -e PASSWORD=pass -v $(pwd):/home/rstudio/sealice rocker/rstudio
```
Istall any required packages using the following command:

```bash
docker exec <container-id> R -e "install.packages('<package-name>')"
```
Where `<package-name>` is the name of the package you want to install and `<container-id>` is the id of the container. The container id can be found by running `docker ps`.

For example to install packages "ggplot2" I would run the following commands:

```bash
docker exec -it $(docker ps -q) R -e 'install.packages(c("ggplot2", "xts", "car", "lmtest", "glmmTMB", "gridExtra"))'
```
Where `docker ps -q` is the container id.

RStudio can be accessed at `http://localhost:8787`
User name is `rstudio` and password is `pass`

## 2. Develop in RStudio
## 3. Containerise the code

Notice that ```docker exec``` is used to install packages in the container. This is not a good practice, as all packages will be lost when the container is stopped and removed. Instead, we should create a Dockerfile and build an image with all the required packages. 
