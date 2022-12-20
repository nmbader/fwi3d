FROM ubuntu:18.04 as builder
MAINTAINER nmbader@sep.stanford.edu
RUN apt-get -y update
RUN apt-get -y install build-essential
RUN apt-get -y install libelf-dev libffi-dev
RUN apt-get -y install pkg-config
RUN apt-get -y install wget git gcc g++ gfortran make cmake vim lsof
RUN apt-get -y install flex libxaw7-dev
RUN apt-get -y install libfftw3-3 libfftw3-dev libssl-dev

RUN apt-get -y  install python3-pip
RUN python3 -m pip install --no-cache-dir --upgrade pip

RUN python3 -m pip install --no-cache-dir numpy &&\
    python3 -m pip install --no-cache-dir jupyter &&\
    python3 -m pip install --no-cache-dir scipy &&\
    python3 -m pip install --no-cache-dir matplotlib

RUN mkdir -p /home
RUN mkdir -p /opt/fwi3d
RUN mkdir -p /home/fwi3d
WORKDIR /home

RUN cd /opt/ &&\
    git clone https://github.com/LLNL/zfp.git &&\
    cd zfp &&\
    mkdir -p build &&\
    cd build &&\
    cmake -DBUILD_SHARED_LIBS=0 .. &&\
    make &&\
    cd /home

ADD src /home/fwi3d/src
ADD external /home/fwi3d/external
ADD examples /home/fwi3d/examples
ADD CMakeLists.txt /home/fwi3d
ADD LICENSE /home/fwi3d
ADD README.md /home/fwi3d

RUN cd /home/fwi3d &&\
    mkdir -p build &&\
    cd external/SEP &&\
    bash ./buildit.sh &&\
    cd ../../build  &&\
    cmake -DCMAKE_INSTALL_PREFIX=/opt/fwi3d/ -DENABLE_ZFP=1 -DCMAKE_PREFIX_PATH="/opt/zfp/build" ../  &&\
    make -j12  &&\
    make install &&\
    cd ../ &&\
    rm -rf build

RUN apt-get -y clean

ENV HOME=/home 
ENV PATH="/opt/fwi3d/bin:${PATH}"
ENV DATAPATH="/tmp/"
RUN echo 'alias python=python3' >> ~/.bashrc