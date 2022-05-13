#include <ogr_spatialref.h>
#include "MathVector.h"
#include "MathMatrix.h"
#include <iostream>
#include <sstream>
#include <complex>
#include <limits>
#include "correction.h"

using namespace std;
using namespace ksg;
using namespace ksg::utils;
using namespace ksg::utils::math;

constexpr int PHI_E = 0;
constexpr int LAMBDA_E = 1;
constexpr int Q = 2;

typedef MathVector<double,3,true> TargetFunctionArg;
typedef MathVector<double,3,true> TargetFunctionResult;
typedef MathMatrix<double,3,3> TargetFunctionDerivative;

constexpr static inline double N(double phi_earth, double a, double eSqr)
{
    return a/sqrt(1-eSqr*pow(sin(phi_earth),2));
}

constexpr static inline double dNdphi(double phi_earth, double a, double eSqr)
{
    return (a*eSqr*sin(phi_earth)*cos(phi_earth))/pow(1-eSqr*pow(sin(phi_earth),2),3/2);
}

static inline TargetFunctionArg createTargetFunctionArg(double phi_earth, double lambda_earth, double q)
{
    TargetFunctionArg ret;
    ret[PHI_E] = phi_earth;
    ret[LAMBDA_E] = lambda_earth;
    ret[Q] = q;
    return ret;
}

static inline TargetFunctionResult targetFunction(const TargetFunctionArg& arg, double a, double eSqr, double h, double phi_s, double lambda_s, double l)
{
    TargetFunctionResult ret;
    double n = N(arg[PHI_E],a,eSqr);
    double cos_phi_e = cos(arg[PHI_E]);
    double cos_phi_s = cos(phi_s);
    ret[0] = (n+h)*cos_phi_e*cos(arg[LAMBDA_E])+cos_phi_s*cos(lambda_s)*arg[Q]-l;
    ret[1] = (n+h)*cos_phi_e*sin(arg[LAMBDA_E])-cos_phi_s*sin(lambda_s)*arg[Q];
    ret[2] = (n*(1-eSqr)+h)*sin(arg[PHI_E])-sin(phi_s)*arg[Q];
    
    return ret;
}

static inline TargetFunctionDerivative targetFunctionDerivative(const TargetFunctionArg& arg, double a, double eSqr, double h, double phi_s, double lambda_s)
{
    TargetFunctionDerivative ret;
    double n = N(arg[PHI_E],a,eSqr);
    double dndphi = dNdphi(arg[PHI_E],a,eSqr);
    double cos_phi_s = cos(phi_s);
    double cos_phi_e = cos(arg[PHI_E]);
    double sin_phi_e = sin(arg[PHI_E]);
    double sin_lambda_e = sin(arg[LAMBDA_E]);
    double cos_lambda_e = cos(arg[LAMBDA_E]);
    
    ret[0][Q] = cos_phi_s*cos(lambda_s);
    ret[1][Q] = -cos_phi_s*sin(lambda_s);
    ret[2][Q] = -sin(phi_s);
    
    ret[0][LAMBDA_E] = -(n+h)*cos_phi_e*sin_lambda_e;
    ret[1][LAMBDA_E] = (n+h)*cos_phi_e*cos_lambda_e;
    ret[2][LAMBDA_E] = 0;
    
    double tmp = (dndphi * cos_phi_e - (n+h)*sin_phi_e);
    
    ret[0][PHI_E] = cos_lambda_e * tmp;
    ret[1][PHI_E] = sin_lambda_e * tmp;
    ret[2][PHI_E] = (1-eSqr)*dndphi*sin_phi_e+(n*(1-eSqr)+h)*cos_phi_e;
    
    return ret;
}

static inline TargetFunctionResult targetFunctionSquared(const TargetFunctionArg& arg, double a, double eSqr, double h, double phi_s, double lambda_s, double l)
{
    auto ret = targetFunction(arg,a,eSqr,h,phi_s,lambda_s,l);
    
    ret.genericInPlaceOperation([](double a)-> double { return a*a;});
    
    return ret;
}

static inline TargetFunctionDerivative targetFunctionSquaredDerivative(const TargetFunctionArg& arg, double a, double eSqr, double h, double phi_s, double lambda_s,double l)
{
    auto value = targetFunction(arg,a,eSqr,h,phi_s,lambda_s,l);
    auto deriv = targetFunctionDerivative(arg,a,eSqr,h,phi_s,lambda_s);
    
    for(int x = 0; x < 3; x++)
        for(int y = 0; y < 3; y++)
        {
            deriv[x][y] = 2 * value[y] * deriv[x][y];
        }
    
    return deriv;
}


class IterationsLimitException : std::exception
{
public:
    IterationsLimitException() :
    exception() {
    }

};

static inline TargetFunctionArg findTargetFunctionRoot(const TargetFunctionArg& start, double a, double eSqr, double h, double phi_s, double lambda_s,double l, double requiredAccuracy, int iterationsLimit)
{
    TargetFunctionArg current = start;
    TargetFunctionResult result = targetFunctionSquared(current,a,eSqr,h,phi_s,lambda_s,l);
    int counter = 0;
    while(result.abs() >= (requiredAccuracy*requiredAccuracy))
    {
        if(counter >= iterationsLimit)
        {
            throw IterationsLimitException();
        }
        auto step = TargetFunctionDerivative::gaussJordanSolve(targetFunctionSquaredDerivative(current,a,eSqr,h,phi_s,lambda_s,l),
                                                               result*(-1.0));
        current += step;
        result = targetFunctionSquared(current,a,eSqr,h,phi_s,lambda_s,l);
        counter++;
    }
    
    return current;
}

constexpr static inline bool isNoData(double v, double noData)
{
    return noData != noData ? (v!=v) : (v==noData);
}

static pair<double,double> calcGeoCoordsFromPix(double xpix, double ypix,double geotransform[])
{
    if(geotransform == NULL) throw runtime_error("calcGeoCordsFromPix: geotranform is NULL");
    double xgeo = geotransform[0] + geotransform[1] * xpix + geotransform[2] * ypix;
    double ygeo = geotransform[3] + geotransform[4] * xpix + geotransform[5] * ypix;
    return pair<double,double>(xgeo,ygeo);
}

static void saveCoordinatesInRaster(GDALDataset* output, int x, int y, double geoX, double geoY) {
    double coords[] = {geoX, geoY};
    int bands[] = {1,2};

    output->RasterIO(GF_Write, x, y, 1, 1, &coords, 1, 1, GDT_Float64, 2, bands,
                        0,0,0,nullptr);
}

void performCorrection(GDALDataset* input, GDALDataset* output, 
                       int inputBand,OGRSpatialReference& geoSrs,
                       double requiredAccuracy, int iterationsLimit)
{
    
    double a = geoSrs.GetSemiMajor();
    double b = geoSrs.GetSemiMinor();
    double central_m = geoSrs.GetProjParm("central_meridian");
    double satelliteHeight = geoSrs.GetProjParm("satellite_height",35785831);
    double eSqr = (a*a - b*b)/(a*a);
    
    auto geogCS = geoSrs.GetAttrNode("GEOGCS");
    char* elips = nullptr;
    geogCS->exportToWkt(&elips);
    
    
    OGRSpatialReference elipsoid;
    elipsoid.importFromWkt((const char**)&elips);
    
    auto tranformation = OGRCreateCoordinateTransformation(&geoSrs,&elipsoid);
    auto reverseTranformation = OGRCreateCoordinateTransformation(&elipsoid,&geoSrs);
    
    int xSize = input->GetRasterXSize();
    int ySize = input->GetRasterYSize();
    
    double geo[6];
    input->GetGeoTransform(geo);
    
    int hasNoData = 0;
    double noData = input->GetRasterBand(inputBand)->GetNoDataValue(&hasNoData);
    
    for(int y=0; y < ySize; y++)
        for(int x=0; x < xSize; x++)
        {
            double h;
            auto band = inputBand;
            input->RasterIO(GF_Read, x, y, 1, 1, &h, 1, 1, GDT_Float64, 1, &band,
                            0,0,0,nullptr);
            auto geoCoords = calcGeoCoordsFromPix(x+0.5,y+0.5,geo);
            
            if(hasNoData && isNoData(h,noData))
            {
                saveCoordinatesInRaster(output,x,y,geoCoords.first, geoCoords.second);
                continue;
            }
            
            
            auto elipsCoords = geoCoords;
            if(!tranformation->Transform(1,&elipsCoords.first,&elipsCoords.second))
            {
                cerr << "Failed to transform pixel " << x <<","<< y << " to elipsoid coords." << endl;
                saveCoordinatesInRaster(output,x,y,numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());
                continue;
            }
            
            elipsCoords.first -= central_m;
            
            elipsCoords.first *= M_PI/180;
            elipsCoords.second *= M_PI/180;
            
            geoCoords.first /= satelliteHeight;
            geoCoords.second /= satelliteHeight;
            
            
            double z = N(elipsCoords.second,a,eSqr)*(1-eSqr)*sin(elipsCoords.second);
            double q = z/sin(geoCoords.second);
            q/=a;
            
            auto startCoords = createTargetFunctionArg(elipsCoords.second,elipsCoords.first,q);
            
            try
            {
                auto finalCoords = findTargetFunctionRoot(startCoords,1,eSqr,h/a,geoCoords.second,geoCoords.first,1+satelliteHeight/a,requiredAccuracy/a, iterationsLimit);
                double newX = finalCoords[LAMBDA_E], newY = finalCoords[PHI_E];
                newX *= 180/M_PI;
                newY *= 180/M_PI;
                
                newX += central_m;

                reverseTranformation->Transform(1,&newX,&newY);

                saveCoordinatesInRaster(output,x,y,newX, newY);
            }
            catch(IterationsLimitException& ex)
            {
                cerr << "Iterations limit for pixel " << x <<","<< y << " exceeded" << endl;
                saveCoordinatesInRaster(output,x,y,numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());
            }
            
        }
    
    CPLFree(tranformation);
    CPLFree(reverseTranformation);
    
}