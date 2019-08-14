# This Dockerfile constructs a docker image that can be used to convert MEG
# data into BIDS format using MNE
#
# Example build:
#   docker build --no-cache --tag winawerlab/meg-bids `pwd`
#

# We are starting from the ubuntu docker...
FROM nben/neuropythy:latest

# Note the Maintainer.
MAINTAINER Noah C. Benson <nben@nyu.edu>

# Stuff to run as root...
USER root
RUN apt-get update \
 && apt-get install -y xvfb libgl1-mesa-glx curl

# && apt-get install -y libglu1-mesa \
# && apt-get install -y libvtk5-dev \
# && apt-get install -y mayavi2

USER $NB_USER

# Install MNE...
RUN mkdir -p mne \
 && cd mne \
 && curl -o ~/mne-env.yml --remote-name https://raw.githubusercontent.com/mne-tools/mne-python/master/environment.yml \
 && conda update -n base conda

RUN conda install vtk
RUN conda env update --file ~/mne-env.yml

RUN pip install xvfbwrapper mne mne-bids

RUN pip install "PyQT5>=5.10"



USER root

COPY mne_bids_export.py /mne_bids_export.py
RUN chown jovyan /mne_bids_export.py && chmod 755 /mne_bids_export.py

# All else can be left as is.
ENTRYPOINT ["python", "/mne_bids_export.py"]


