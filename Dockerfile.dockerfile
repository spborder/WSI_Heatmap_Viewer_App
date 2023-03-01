# Dockerfile for FUSION: Functional Unit State Identification and Navigation with WSI
FROM python:3.8

LABEL maintainer="Sam Border CMI Lab <samuel.border@medicine.ufl.edu>"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    openslide-tools \
    python3-openslide

ENV RUNTYPE=='AWS'

COPY . ./
RUN python3 -m pip install --upgrade pip

# Terracotta setup
WORKDIR /assets/slide_info/
RUN git clone https://github.com/DHI-GRAS/terracotta.git
RUN python3 -m pip install -r ./requirements.txt --no-cache-dir
EXPOSE 5000

RUN terracotta serve -r full_optimized_tifs/{slide}_{band}.tif --port 5000

# Running the app
WORKDIR /
RUN python3 -m pip install -r ./requirements.txt --no-cache-dir
RUN python3 -m pip freeze

EXPOSE 8000

ENTRYPOINT [ "python3" ]
CMD ["WSI_Heatmap_Viewer.py"]
