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
#include "correction_utils.h"

class GEOSHeightCorrector
{
private:
    double a, b, central_m, satelliteHeight;
    double eSqr;

    OGRSpatialReference geosSrs, elipsoidSrs;

    OGRCoordinateTransformation *transformation, *reverseTranformation;

public:
    GEOSHeightCorrector(OGRSpatialReference &geoSrs);

    ~GEOSHeightCorrector();

    void calculateNewCoordinatesForRaster(
        GDALDataset *input,
        GDALDataset *output,
        int inputBand,
        double requiredAccuracy,
        int iterationsLimit,
        ksg::NumericMethod numericMethod = ksg::NumericMethod::NETWON,
        bool useQuadraticForm = false
    );

    std::pair<double, double> calculateNewCoordinates(
        std::pair<double, double> geoCoords,
        std::pair<double, double> ellipsCoords,
        double objectHeight,
        double requiredAccuracy,
        int iterationsLimit,
        ksg::NumericMethod numericMethod = ksg::NumericMethod::NETWON,
        bool useQuadraticForm = false
    );

    std::pair<double, double> transformToEllipsCoordinates(
        std::pair<double, double> geoCoords);
    std::pair<double, double> transformToGeosCoordinates(
        std::pair<double, double> ellipsCoords);
};

#endif /* CORRECTION_H */
