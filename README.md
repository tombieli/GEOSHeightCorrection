# GEOS Height Correction

## Introduction

This software is a supplement to the publication "A Parallax Shift Effect Correction Based on Cloud Height for Geostationary Satellites and Radar Observations" *[[1]](#references)*.

This program performs correction of cloud's parallax shift using its a priori known height. Program is intended to use with raster images in [Geostationary projection](https://proj.org/operations/projections/geos.html) *[[2]](#references)*.

### C++ program

Repository contains C++ program `geosheightcorrection` which calculates new coordinates  in geostationary view of every pixel using its original coordinates (geostationary coordinates of cloud's image) and a priori known height of cloud. Usually cloud height can be calculated using cloud top temperature or in more sophisticated way - like CTH product provided by Eumetcast *[[3]](#references)*.

C++ program requires following libraries to run:
* GDAL *[[4]](#references)*,
* boost-program_options *[[5]](#references)*

### Wrapping script

This repository also contains script `perform_geos_height_correction.sh` which allows to correct existing images in Geostationary Projection using priori known height of clouds. Script produces raster in which data connected with cloud top are moved to relevant nadir location.

To run, wrapping script requires above C++ program present in `PATH` and installed GDAL *[[4]](#references)* utilities including python scripts.

## Build

Project can be build using Docker or CMake.

### Docker

To build project on docker please use following command.

```shell
docker build -t geosheightcorrection .
```

Build docker image will contain both `geosheightcorrection` program and `perform_geos_height_correction.sh` script in `PATH`. For detailed information about usage please read following [section](#usage).

### CMake (Linux)

To build C++ program on Linux it following elements are required:
* C++ compiler with support of C++11 a least (g++ is suggested),
* CMake program,
* make program (or other tool supported by CMake),
* GDAL library with headers,
* Boost program_options library with headers

To build C++ program with CMake please execute following commands in main repository directory:

```shell
cmake build
cd build
make
```

After those steps in `build` directory should reside program `geosheightcorrection` ready to use.

## Usage

### Assumptions

This correction software can be used either from docker or directly by operation system.

#### Docker

For docker usage it is assumed that command is preceded by docker prefix:

Powershell on Windows

```powershell
docker run --rm -v "${pwd}/data:/app/data" -ti geosheightcorrection <executable_name> <executable_parameters>
```

On Linux:

```shell
docker run --rm -v `pwd`/data:/app/data -ti geosheightcorrection <executable_name> <executable_parameters>
```

`pwd` should point to main repository directory.

#### Linux

In case of running correction executable directly in operating system it is assumed that GDAL utilities and `geosheightcorrection` are present in `PATH`.

### Sample Data

`data` directory contains sample data that could be used for evaluation of C++ program and wrapping script.

Those data is following:
* `ctt_201507251300.tif` - cloud top temperature extracted from SEVIRI image (channel IR 10.8) [[6]](#references),
* `h_201507251300.tif` - cloud top height calculated with simple algorithm using data from `ctt_201507251300.tif`.

### C++ program

C++ program run without parameters:

```shell
geosheightcorrection
```

will display option description:

```text
GEOS Height correction options:
  --help                        Shows help message
  --input arg                   Location of image containing height information.                     
  --output arg                  Location of result raster with corrected coordiantes
  --height-band arg (=1)        Band with height information [m]. 
  --requierd-accuracy arg (=10) Required accuracy [m].
  --iterations-limit arg (=100) Maximum number of iteration of netwon method per pixel.
```

To calculate corrected GEOS coordinates using raster `h_201507251300.tif` please execute following command:

```shell
geosheightcorrection --input data/h_201507251300.tif --output data/new_coordinates.tif
```

After that operation in `data` directory new file `new_coordinates.tif` should be present. This raster will contain two `FLOAT_64` bands. Those band will represent new coordinates (x and y respectively) in geostationary source reference system (same as for input and output image).

### Wrapping script

Wrapping script can perform more complex processing which produces fully adjusted image with approximately same content.

Wrapping script has following options:
* `--input` - path to image requiring correction,
* `--height` - path to image containing cloud height data (in meters),
* `--output` - location where corrected image will be saved,
* `--height-band` - number of band containing height information in height raster. If not specified `1` is assumed,
* `--algorithms` - list of interpolation algorithms for each band respectively. If not defined `average` is assumed for each band. For bands which contain quantified data by definition `nearest` algorithm is suggested. For more information please refer to help of `gdal_grid` program.

To run full correction and resampling of sample data please run following command:

```shell
perform_geos_height_correction.sh --input data/ctt_201507251300.tif --height data/h_201507251300.tif --output data/ctt_201507251300_corrected.tif 
```

After that command corrected raster `ctt_201507251300_corrected.tif` should appear in `data` directory.

Alternatively raster can be corrected with different interpolation algorithm

```shell
perform_geos_height_correction.sh --input data/ctt_201507251300.tif --height data/h_201507251300.tif --output data/ctt_201507251300_corrected_nearest.tif --algorithms "nearest"
```

## References

* [1] [Bieliński, T. (2020). A Parallax Shift Effect Correction Based on Cloud Height for Geostationary Satellites and Radar Observations. Remote Sensing, 12, 1-20. https://doi.org/10.3390/rs12030365](https://doi.org/10.3390/rs12030365)
* [2] [PROJ contributors (2020). PROJ coordinate transformation software library. Open Source Geospatial Foundation. URL https://proj.org/. DOI: 10.5281/zenodo.5884394](https://proj.org/)
* [3] [EUMETSAT. MSG Meteorological Products Extraction Facility
Algorithm Specification Document. (Access 2022-05-20)](https://www.eumetsat.int/media/38993)
* [4] [GDAL/OGR contributors (2022). GDAL/OGR Geospatial Data Abstraction software Library. Open Source Geospatial Foundation. URL https://gdal.org DOI: 10.5281/zenodo.5884351](https://gdal.org)
* [5] [Boost program options library. (Access 2022-05-20)](https://www.boost.org/doc/libs/1_79_0/doc/html/program_options.html)
* [6] [MSG 0 degree description on EUMETSAT page. (Access 2022-05-20)](https://eumetview.eumetsat.int/static-images/MSG/)
