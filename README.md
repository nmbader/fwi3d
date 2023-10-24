# FWI3D: **_C++_** library for 3D modeling and inversion of seismic waves

- [Description](#Description)
- [Prerequisites](#Prerequisites)
- [Installation](#Installation)
- [Data format](#Data-format)

## Description

This library performs three-dimensional modeling and inversion of seismic data using fourth-order summation-by-parts finite-difference operators on a regular cartesian grid. Three main configurations are available: acoustic (variable density), elastic isotropic, elastic VTI.

Key features include:

* incorporated linear and non-linear solvers (including non-linear CG and l-BFGS)
* many practical options for full-waveform inversion (FWI)
* **_MPI_** support for parallelization over shots

## Installation using Docker

Start by cloning the current repository
```
git clone https://github.com/nmbader/fwi3d.git
cd fwi3d
```

Build the docker image (it should take a few minutes)
```
docker build -f Dockerfile -t fwi3d .
```

Run a container
```
docker run -it -p 8080:8080 fwi3d
```

By default a bash shell will be opened at /home inside the container.
Run jupyter notebook from within the container
```
jupyter notebook --ip 0.0.0.0 --port 8080 --no-browser --allow-root &
```

Open the browser at *localhost:8080/​* and use the printed token above to authenticate.

The image will build the **fwi3d** library in single precision.

## Installation without a Docker

### Prerequisites

The **fwi3d** library is written entirely in **_C++_** and has been built on Linux environments (**centos 7**) using **_CMAKE_**. The **_FFTW3_** library is required. Moreover, **_Python_** and **_Jupyter_** are only required for data format conversion and examples generation.

### Installation

Clone and install the **fwi3d** library
```
# get the code
git clone https://github.com/nmbader/fwi3d.git
cd fwi3d

# create subdirectories
mkdir build local

# Build the external SEPlib library needed for the IO
cd external/SEP
bash ./buildit.sh

# When it is done, build the main library
cd ../../build
cmake -DCMAKE_INSTALL_PREFIX=../local ../
make -j12
make install

# clean up the build directory
rm -rf *
```
By default, the **fwi3d** library is built in single precision. For double precision, add the flag `-DENABLE_DOUBLE_PRECISION=1`
 to the **_cmake_** command. 

All executables (and python scripts) will be installed into the subdirectory *fwi3d/local/bin*. It will be more convenient to add this path to the environment variable `PATH` in order to run the examples seeminglessly.

**_MPI_** is also available for parallelization over seismic sources. If **_CMAKE_** cannot locate an **_MPI_** installation automatically, add the corresponding path manually by setting the flag `-DCMAKE_PREFIX_PATH=path_to_mpi_directory`. Otherwise, **_MPI_** will be deactivated.


## Data format

By default, the main executables read and write data in *SEPlib* format (with little-endian binaries). Alternatively, native binary format is also accepted provided that a description file is built. A python script is provided to convert to/from *SEPlib* from/to **_numpy_**, and a python class to write and read directly from python to *SEPlib*.

Refer to the [examples](https://github.com/nmbader/fwi3d/tree/master/examples) for a simple modeling test. Refer to the 2D version of this software [examples](https://github.com/nmbader/fwi2d) with more modeling and inversion examples.
