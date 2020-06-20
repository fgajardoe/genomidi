FROM debian:testing
LABEL maintainer Felipe Gajardo <fgajardoe@gmail.com>

# Install software.
RUN apt-get update && \
    apt-get install -y \
     gcc \
     g++ \
     git \
     cmake \
     timidity \
     bedtools \
     vim \
     r-base \
     libcurl4-openssl-dev \
     fluid-soundfont-gm \
    && \
    Rscript -e 'install.packages(c("BiocManager"))' && \
    Rscript -e 'BiocManager::install("GenomicRanges")' && \
    cd /tmp && \
    git clone https://github.com/markc/midicomp && \
    mkdir /tmp/midicomp/build && \
    cd /tmp/midicomp/build && \
    cmake .. && \
    make && \
    mv /tmp/midicomp/build/midicomp /usr/bin/. 

# Add /opt/bin to PATH.
ENV PATH /opt/bin:$PATH

# Set user.
RUN useradd -ms /bin/bash genomidi
RUN echo 'genomidi:genomidi' | chpasswd
USER genomidi
WORKDIR /home/genomidi

RUN alias ll='ls -l'
RUN alias l='ls'
RUN alias xx='exit'
RUN mkdir /home/genomidi/data

ADD symphony.R /home/genomidi/symphony.R
ADD genomidi.R /home/genomidi/genomidi.R


CMD ["/bin/bash"]
