#include "fixtures.h"
#include <boost/test/unit_test.hpp>

namespace ksg
{

    const char* Fixture::proj = "+proj=geos +lon_0=0 +h=35785831 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs";

    OGRSpatialReference Fixture::initSrs()
    {
        OGRSpatialReference srs;
        BOOST_TEST(srs.importFromProj4(proj) == CE_None);
        return srs;
    }

    Fixture::Fixture()
        : srs(initSrs()), corrector(srs)
    {
        a = srs.GetSemiMajor();
        b = srs.GetSemiMinor();
        central_m = srs.GetProjParm("central_meridian");
        satelliteHeight = srs.GetProjParm("satellite_height", 35785831);
        eSqr = (a * a - b * b) / (a * a);
    }

    PixelFixture::PixelFixture()
        : Fixture(), geotransform({5.5702484773397446e+06, -3.0004031658172607e+03, 0.0000000000000000e+00, -5.5702484773397446e+06, 0.0000000000000000e+00, 3.0004031658172607e+03})
    {
    }

}