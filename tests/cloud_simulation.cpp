#include "cloud_simulation.h"
#include <cmath>
#include "geodetic_utils.h"

using namespace std;

namespace ksg
{

    Point3D calculateCloudPosition(double lat, double lon, double height, double a, double eSqr, double lambda0)
    {
        Point3D ret;
        auto n = N(lat, a, eSqr);
        auto lonRed = lon - lambda0;
        ret(X) = (n + height) * cos(lat) * cos(lonRed);
        ret(Y) = (n + height) * cos(lat) * sin(lonRed);
        ret(Z) = (n * (1 - eSqr) + height) * sin(lat);

        return ret;
    }

    std::pair<double, double> calculateGEOSCorrdsFromXYZ(const Point3D &point, double a, double h)
    {
        double l = a + h;
        double phi = atan(point(Z) / (sqrt(pow(point(X) - l, 2) + pow(point(Y), 2))));
        double lambda = -atan(point(Y) / (point(X) - l));

        return make_pair(lambda * h, phi * h);
    }
}