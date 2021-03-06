/* 
 * File:   correction.h
 * Author: tombieli
 *
 * Created on 1 August 2016, 15:47
 */

#ifndef CORRECTION_H
#define CORRECTION_H

#include <ogr_spatialref.h>
#include <gdal_priv.h>

void performCorrection(GDALDataset* input, GDALDataset* output,
                       int inputBand,OGRSpatialReference& geoSrs,double requiredAccuracy, int iterationsLimit);

#endif /* CORRECTION_H */

