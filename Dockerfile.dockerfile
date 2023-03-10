# Dockerfile for FUSION: Functional Unit State Identification and Navigation with WSI
FROM python:3.6

LABEL maintainer="Sam Border CMI Lab <samuel.border@medicine.ufl.edu>"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    openslide-tools \
    python3-openslide

ENV RUNTYPE='AWS'
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal

COPY . ./
RUN python3 -m pip install --upgrade pip

# Running the app
WORKDIR /
RUN python3 -m pip install -r ./requirements.txt --no-cache-dir
RUN python3 -m pip freeze > pip_installed_packages.txt

EXPOSE 8000

ENTRYPOINT [ "python3" ]
CMD ["WSI_Heatmap_Viewer.py"]
