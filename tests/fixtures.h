#include <correction.h>
#include <georeference_utils.h>
#include <ogr_spatialref.h>

namespace ksg
{

    struct Fixture
    {
        static const char *proj;

        OGRSpatialReference srs;
        GEOSHeightCorrector corrector;
        double a, b, central_m, satelliteHeight;
        double eSqr;

        static OGRSpatialReference initSrs();
        Fixture();
    };

    struct PixelFixture : public Fixture
    {
        ksg::Geotransform<> geotransform;

        PixelFixture();
    };
}