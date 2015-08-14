Multi-file zonal statistics Python / FME modules
------------------------------------------------

This repository contains code for calculating zonal statistics efficiently across a large number of geographically matched input rasters and a single zones dataset. 

The python code is derived from the python [rasterstats](https://github.com/perrygeo/python-rasterstats) module. This has been modified to accept a list of input raster files and to store the rasterized zones data in memory and apply it to each data raster in turn without the need to re-rasterize.

Although it is not necessary to use it this way, the code is developed for use via FME through a custom transformer. This custom transformer has been provided, along with a trivial demonstration workbench to illustrate its use. The envisaged use for this within MAP is to produce and summarise data tables from model output data cubes by using FME to process the raster stats outputs as they are produced. 

Note that you will either need to install the requisite python modules (gdal etc) into the copy of Python that FME is using, or else set FME to use a different python interpreter on your machine that already has them installed.

