# Dockerfile for Terracotta tile server
FROM continuumio/miniconda3

LABEL maintainer = "Sam Border CMI Lab <samuel.border@medicine.ufl.edu>" 
WORKDIR /

# Getting bash commands
SHELL ["/bin/bash", "--login", "-c"]

# Copying necessary files
COPY environment.yml .
COPY terracotta .
COPY full_optimized_tifs .

# Creating environment
RUN conda env create -f environment.yml

# Initializing conda in bash configs
RUN conda init bash
# Activating environment
RUN conda activate terracotta

# Exposing port for tile-server
EXPOSE 5000

# Entrypoint commands
ENTRYPOINT ["terracotta serve -r full_optimized_tifs/{slide}_{band}.tif --port 5000"]