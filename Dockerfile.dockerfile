# Dockerfile for WSI_Heatmap_Viewer Dash app
# Not sure if it will work :/

FROM python:3.6

LABEL maintainer="Sam Border CMI Lab <samuel.border@medicine.ufl.edu>"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    openslide-tools \
    python3-openslide

RUN git clone https://github.com/spborder/WSI_Heatmap_Viewer_App.git

COPY . ./
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install -r ./requirements.txt --no-cache-dir
RUN python3 -m pip install tiffslide --no-cache-dir
RUN python3 -m pip freeze


CMD gunicorn -b 0.0.0:80 app.app:server
# docker run -p {port#}:80 {tag-name}