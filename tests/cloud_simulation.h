#include <dlib/geometry/vector.h>
#include <utility>

namespace ksg {

    typedef dlib::vector<double, 3> Point3D;
    constexpr unsigned long X = 0;
    constexpr unsigned long Y = 1;
    constexpr unsigned long Z = 2;

    Point3D calculateCloudPosition(double lat, double lon, double height, double a, double eSqr, double lambda0);
    std::pair<double, double> calculateGEOSCorrdsFromXYZ(const Point3D& point, double a, double h);

}