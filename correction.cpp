#include <ogr_spatialref.h>
#include <iostream>
#include <sstream>
#include <complex>
#include <limits>
#include <functional>
#include <dlib/geometry/vector.h> 
#include <dlib/matrix.h> 
#include "correction.h"
#include "georeference_utils.h"
#include "correction_utils.h"
#include "geodetic_utils.h"


using namespace std;
using namespace std::placeholders;
using namespace ksg;
using namespace dlib;

constexpr int PHI_E = 0;
constexpr int LAMBDA_E = 1;
constexpr int Q = 2;

typedef dlib::vector<double, 3> TargetFunctionArg;
typedef dlib::vector<double, 3> TargetFunctionResult;
typedef matrix<double, 3, 3> TargetFunctionDerivative;
typedef matrix<double, 3, 3> TargetFunctionHessian;


static inline TargetFunctionArg createTargetFunctionArg(double phi_earth, double lambda_earth, double q)
{
    TargetFunctionArg ret;
    ret(PHI_E) = phi_earth;
    ret(LAMBDA_E) = lambda_earth;
    ret(Q) = q;
    return ret;
}


class ParallaxProblem {

    double a;
    double eSqr;
    double h;
    double phi_s; 
    double lambda_s; 
    double l;

public:
    ///
    /// @param a Ellipsoid semi-major axis lenght
    /// @param eSqr Ellipsoid flattering coefficient squared
    /// @param h Object height above surface
    /// @param phi_s Satellite vertical view angle in radians
    /// @param lambda_s Satellite horizontal view angle in radians
    /// @param l Satellite distance from center of ellipsoid
    ParallaxProblem(double a, double eSqr, double h, double phi_s, double lambda_s, double l)
    : a(a), eSqr(eSqr), h(h), phi_s(phi_s), lambda_s(lambda_s), l(l)
    {}

    inline double targetFunctionScalarValue(const TargetFunctionArg &arg) const {
        return targetFunctionVectorValue(arg).length();
    }

    inline TargetFunctionResult targetFunctionVectorValue(const TargetFunctionArg &arg) const
    {
        TargetFunctionResult ret;
        double n = N(arg(PHI_E), a, eSqr);
        double cos_phi_e = cos(arg(PHI_E));
        double cos_phi_s = cos(phi_s);
        ret(0) = (n + h) * cos_phi_e * cos(arg(LAMBDA_E)) + cos_phi_s * cos(lambda_s) * arg(Q) - l;
        ret(1) = (n + h) * cos_phi_e * sin(arg(LAMBDA_E)) - cos_phi_s * sin(lambda_s) * arg(Q);
        ret(2) = (n * (1 - eSqr) + h) * sin(arg(PHI_E)) - sin(phi_s) * arg(Q);

        return ret;
    }

    inline TargetFunctionDerivative targetFunctionJacobian(const TargetFunctionArg &arg) const
    {
        TargetFunctionDerivative ret;
        double n = N(arg(PHI_E), a, eSqr);
        double dndphi = dNdphi(arg(PHI_E), a, eSqr);
        double cos_phi_s = cos(phi_s);
        double cos_phi_e = cos(arg(PHI_E));
        double sin_phi_e = sin(arg(PHI_E));
        double sin_lambda_e = sin(arg(LAMBDA_E));
        double cos_lambda_e = cos(arg(LAMBDA_E));

        ret(0, Q) = cos_phi_s * cos(lambda_s);
        ret(1, Q) = -cos_phi_s * sin(lambda_s);
        ret(2, Q) = -sin(phi_s);

        ret(0, LAMBDA_E) = -(n + h) * cos_phi_e * sin_lambda_e;
        ret(1, LAMBDA_E) = (n + h) * cos_phi_e * cos_lambda_e;
        ret(2, LAMBDA_E) = 0;

        double tmp = (dndphi * cos_phi_e - (n + h) * sin_phi_e);

        ret(0, PHI_E) = cos_lambda_e * tmp;
        ret(1, PHI_E) = sin_lambda_e * tmp;
        ret(2, PHI_E) = (1 - eSqr) * dndphi * sin_phi_e + (n * (1 - eSqr) + h) * cos_phi_e;

        return ret;
    }

    inline TargetFunctionHessian targetFunctionHessian(const TargetFunctionArg &arg) const {
        auto jacobi = targetFunctionJacobian(arg);
        return trans(jacobi)*jacobi;
    }

    ///
    /// @param value target function value
    /// @param accuracy required accuracy
    /// @return true if provided target function value satisfies required accuracy
    static inline bool doesSatsifyRequiredAccuracy(const TargetFunctionResult value, double accuracy) {
        return value.length() <= accuracy;
    }

};

class ParallaxProblemSquared {
    ParallaxProblem original;

public:
    ///
    /// @param a Ellipsoid semi-major axis lenght
    /// @param eSqr Ellipsoid flattering coefficient squared
    /// @param h Object height above surface
    /// @param phi_s Satellite vertical view angle in radians
    /// @param lambda_s Satellite horizontal view angle in radians
    /// @param l Satellite distance from center of ellipsoid
    ParallaxProblemSquared(double a, double eSqr, double h, double phi_s, double lambda_s, double l)
    : original(a,eSqr, h, phi_s, lambda_s, l)
    {}

    inline double targetFunctionScalarValue(const TargetFunctionArg &arg) const {
        return targetFunctionVectorValue(arg).length();
    }

    inline TargetFunctionResult targetFunctionVectorValue(const TargetFunctionArg &arg) const
    {
        auto ret = original.targetFunctionVectorValue(arg);

        for(auto& elem : ret) {
            elem = elem * elem;
        } 

        return ret;
    }

    inline TargetFunctionDerivative targetFunctionJacobian(const TargetFunctionArg &arg) const
    {
        auto value = original.targetFunctionVectorValue(arg);
        auto deriv = original.targetFunctionJacobian(arg);

        for (int fun = 0; fun < 3; fun++)
            for (int var = 0; var < 3; var++)
            {
                                      //    v Should be var, but does not work
                deriv(fun, var) = 2 * value(fun) * deriv(fun, var);
            }

        return deriv;
    }

    inline TargetFunctionHessian targetFunctionHessian(const TargetFunctionArg &arg) const {
        auto jacobi = targetFunctionJacobian(arg);
        return trans(jacobi)*jacobi;
    }

    ///
    /// @param value target function value
    /// @param accuracy required accuracy
    /// @return true if provided target function value satisfies required accuracy
    static inline bool doesSatsifyRequiredAccuracy(const TargetFunctionResult value, double accuracy) {
        return value.length() <= pow(accuracy, 2)/2;
    }
};


class IterationsLimitException : std::exception
{
public:
    IterationsLimitException() : exception()
    {
    }
};


template<class P>
static inline TargetFunctionArg findTargetFunctionRootByLM(const TargetFunctionArg &start, const P& problem, double requiredAccuracy, int iterationsLimit)
{
    TargetFunctionArg current = start;
    auto result = problem.targetFunctionVectorValue(current);
    auto hessian = problem.targetFunctionHessian(current);
    int counter = 0;
    double damping = max(hessian);
    const auto I = identity_matrix<double, 3>();
    while (!P::doesSatsifyRequiredAccuracy(result, requiredAccuracy))
    {
        auto jacobi = problem.targetFunctionJacobian(current);
        TargetFunctionArg grad = (trans(jacobi) * result);
        hessian = problem.targetFunctionHessian(current);
        bool stepSucceed = false;
        do {
            TargetFunctionArg step = inv(hessian + damping * I) * grad;
            auto next = current - step;
            auto nextResult = problem.targetFunctionVectorValue(next);
            stepSucceed = nextResult.length() < result.length();
            if(!stepSucceed) {
                damping *= 2;
            }
            else {
                damping /= 2;
                current = next;
                result = nextResult;
            }
            counter++;
            if (counter >= iterationsLimit)
            {
                throw IterationsLimitException();
            }
        }
        while(!stepSucceed);
        
    }

    return current;
}

template<class P>
static inline TargetFunctionArg findTargetFunctionRootByNewton(const TargetFunctionArg &start, const P& problem, double requiredAccuracy, int iterationsLimit)
{
    TargetFunctionArg current = start;
    auto result = problem.targetFunctionVectorValue(current);
    int counter = 0;
    double alpha = 1;
    while (!P::doesSatsifyRequiredAccuracy(result, requiredAccuracy))
    {
        auto jacobi = problem.targetFunctionJacobian(current);
        auto invJacobi = inv(jacobi);
        TargetFunctionArg step = invJacobi * result;
        bool stepSucceed = false;
        do {
            auto next = current - alpha * step;
            auto nextResult = problem.targetFunctionVectorValue(next);
            stepSucceed = nextResult.length() < result.length();
            if(!stepSucceed) {
                alpha /= 2;
            }
            else {
                alpha *= 2;
                current = next;
                result = nextResult;
            }
            counter++;
            if (counter >= iterationsLimit)
            {
                throw IterationsLimitException();
            }
        }
        while(!stepSucceed);
        
    }

    return current;
}

using parallaxSolver = 
    std::function<TargetFunctionArg (
                            const TargetFunctionArg &start,
                            double a, double eSqr, double h, double phi_s, double lambda_s, double l,
                            double requiredAccuracy, int iterationsLimit
    )>;

template<class P>
using numericMethod = 
    std::function<TargetFunctionArg (const TargetFunctionArg &start, const P& problem, double requiredAccuracy, int iterationsLimit)>;


template<class P>
static inline numericMethod<P> getNumericMethodImplememntation(NumericMethod method) {
    switch(method) {
        case NumericMethod::LEVENBERG_MARQUARD:
            return findTargetFunctionRootByLM<P>;
        case NumericMethod::NETWON:
            return findTargetFunctionRootByNewton<P>;
        default:
            stringstream str;
            str << method;
            throw std::logic_error("Unknown numeric method: "+str.str());
    }
}


static inline parallaxSolver prepareParallaxSolver(NumericMethod method, bool squared) {
    if(!squared) {
        return [=](const TargetFunctionArg &start,
                    double a, double eSqr, double h, double phi_s, double lambda_s, double l,
                    double requiredAccuracy, int iterationsLimit) {
                        ParallaxProblem p{a,eSqr,h,phi_s, lambda_s, l};
                        auto nm = getNumericMethodImplememntation<ParallaxProblem>(method);
                        return nm(start,p, requiredAccuracy, iterationsLimit);
                    };
    }
    else {
        return [=](const TargetFunctionArg &start,
            double a, double eSqr, double h, double phi_s, double lambda_s, double l,
            double requiredAccuracy, int iterationsLimit) {
                ParallaxProblemSquared p{a,eSqr,h,phi_s, lambda_s, l};
                auto nm = getNumericMethodImplememntation<ParallaxProblemSquared>(method);
                return nm(start,p, requiredAccuracy, iterationsLimit);
            };
    }
}

constexpr static inline bool isNoData(double v, double noData)
{
    return noData != noData ? (v != v) : (v == noData);
}

static void saveCoordinatesInRaster(GDALDataset *output, int x, int y, double geoX, double geoY)
{
    double coords[] = {geoX, geoY};
    int bands[] = {1, 2};

    if (output->RasterIO(GF_Write, x, y, 1, 1, &coords, 1, 1, GDT_Float64, 2, bands,
                         0, 0, 0, nullptr) != CE_None)
    {
        cerr << "Error while writing pixel " << x << ", " << y << " to output raster." << endl;
    }
}

GEOSHeightCorrector::GEOSHeightCorrector(OGRSpatialReference &geosSrs)
    : geosSrs(geosSrs)
{
    a = geosSrs.GetSemiMajor();
    b = geosSrs.GetSemiMinor();
    central_m = geosSrs.GetProjParm("central_meridian");
    satelliteHeight = geosSrs.GetProjParm("satellite_height", 35785831);
    eSqr = (a * a - b * b) / (a * a);

    auto geogCS = geosSrs.GetAttrNode("GEOGCS");
    char *elips = nullptr;
    geogCS->exportToWkt(&elips);

    elipsoidSrs.importFromWkt(elips);
    CPLFree(elips);

    transformation = OGRCreateCoordinateTransformation(&geosSrs, &elipsoidSrs);
    reverseTranformation = OGRCreateCoordinateTransformation(&elipsoidSrs, &geosSrs);
}

GEOSHeightCorrector::~GEOSHeightCorrector()
{
    OGRCoordinateTransformation::DestroyCT(transformation);
    OGRCoordinateTransformation::DestroyCT(reverseTranformation);
}

std::pair<double, double> GEOSHeightCorrector::calculateNewCoordinates(
    std::pair<double, double> geoCoords,
    std::pair<double, double> elipsCoords,
    double objectHeight,
    double requiredAccuracy,
    int iterationsLimit,
    ksg::NumericMethod numericMethod,
    bool useQuadraticForm
)
{
    elipsCoords.first -= central_m;

    elipsCoords.first *= M_PI / 180;
    elipsCoords.second *= M_PI / 180;

    geoCoords.first /= satelliteHeight;
    geoCoords.second /= satelliteHeight;

    double distanceFormSatellite;

    if(geoCoords.second > numeric_limits<double>::epsilon()) {
        double z = N(elipsCoords.second, a, eSqr) * (1 - eSqr) * sin(elipsCoords.second);
        distanceFormSatellite = z / sin(geoCoords.second);
    } else {
        distanceFormSatellite = satelliteHeight;
    }


    double q = distanceFormSatellite - 2*objectHeight;
    q /= a;

    auto startCoords = createTargetFunctionArg(elipsCoords.second, elipsCoords.first, q);

    try
    {
        auto solver = prepareParallaxSolver(numericMethod, useQuadraticForm);
        auto finalCoords = solver(startCoords, 1, eSqr, objectHeight / a, geoCoords.second, geoCoords.first, 1 + satelliteHeight / a, requiredAccuracy / a, iterationsLimit);
        double newX = finalCoords(LAMBDA_E), newY = finalCoords(PHI_E);
        newX *= 180 / M_PI;
        newY *= 180 / M_PI;

        newX += central_m;

        return std::make_pair(newX, newY);
    }
    catch (IterationsLimitException &ex)
    {
        cerr << "Iterations limit for geos " << geoCoords.first << "," << geoCoords.second << " exceeded" << endl;
        return std::make_pair(numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());
    }
}

std::pair<double,double> GEOSHeightCorrector::transformToEllipsCoordinates(
    std::pair<double,double> geoCoords
) {
    auto elipsCoords = geoCoords;
    if (!transformation->Transform(1, &elipsCoords.first, &elipsCoords.second))
    {
        cerr << "Failed to transform " << geoCoords.first << ", " << geoCoords.second << " from geos to elipsoid coords.\n";
        return std::make_pair(numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());
    }

    // if(!elipsoidSrs.EPSGTreatsAsLatLong()) {
    //     swap(elipsCoords.first, elipsCoords.second);
    // }

    return elipsCoords;
}

std::pair<double,double> GEOSHeightCorrector::transformToGeosCoordinates(
    std::pair<double,double> elipsCoords
) {
    // if(!elipsoidSrs.EPSGTreatsAsLatLong()) {
    //     swap(elipsCoords.first, elipsCoords.second);
    // }
    auto geoCoords = elipsCoords;
    if (!reverseTranformation->Transform(1, &geoCoords.first, &geoCoords.second))
    {
        cerr << "Failed to transform " << elipsCoords.first << ", " << elipsCoords.second << " from elipsoid to geos coords.\n";
        return std::make_pair(numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());
    }

    return geoCoords;
}


void GEOSHeightCorrector::calculateNewCoordinatesForRaster(
    GDALDataset *input, GDALDataset *output,
    int inputBand, double requiredAccuracy, int iterationsLimit,
    ksg::NumericMethod numericMethod,
    bool useQuadraticForm
)
{
    int xSize = input->GetRasterXSize();
    int ySize = input->GetRasterYSize();

    Geotransform<> geo;
    input->GetGeoTransform(geo.geotransform);

    int hasNoData = 0;
    double noData = input->GetRasterBand(inputBand)->GetNoDataValue(&hasNoData);

    for (int y = 0; y < ySize; y++)
        for (int x = 0; x < xSize; x++)
        {
            double h;
            auto band = inputBand;
            if (input->RasterIO(GF_Read, x, y, 1, 1, &h, 1, 1, GDT_Float64, 1, &band,
                                0, 0, 0, nullptr) != CE_None)
            {
                cerr << "Failed to fetch pixel " << x << "," << y << " from input raster." << endl;
                saveCoordinatesInRaster(output, x, y, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());
                continue;
            }
            auto geoCoords = geo.calcGeoCoordsFromPix(x + 0.5, y + 0.5);

            if (hasNoData && isNoData(h, noData))
            {
                saveCoordinatesInRaster(output, x, y, geoCoords.first, geoCoords.second);
                continue;
            }

            auto elipsCoords = geoCoords;
            if (!transformation->Transform(1, &elipsCoords.first, &elipsCoords.second))
            {
                cerr << "Failed to transform pixel " << x << "," << y << " from geos to elipsoid coords." << endl;
                saveCoordinatesInRaster(output, x, y, numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());
                continue;
            }

            auto result = calculateNewCoordinates(geoCoords, elipsCoords, h, requiredAccuracy, iterationsLimit, numericMethod, useQuadraticForm);

            if (isnan(result.first) || isnan(result.second))
            {
                cerr << "Failed to correct pixel" << x << "," << y << "." << endl;
            }
            else
            {
                reverseTranformation->Transform(1, &result.first, &result.second);
            }
            saveCoordinatesInRaster(output, x, y, result.first, result.second);
        }
}