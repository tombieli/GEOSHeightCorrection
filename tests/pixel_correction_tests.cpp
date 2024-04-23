#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <vector>
#include <cmath>

#include "fixtures.h"
#include "cloud_simulation.h"
#include "test_utils.h"

using namespace std;
using namespace ksg;
namespace data = boost::unit_test::data;

static vector<PixPosition> positions = {
    {1856,1856}, //Center
    {1832, 51}, //Suspicious South East
    {1488, 3439} //Gdańsk
  };

static vector<int> iterations = {100, 1000};
static vector<double> requiredAccuracies = {10, 1};
static vector<NumericMethod> numericMethods = {NumericMethod::NETWON, NumericMethod::LEVENBERG_MARQUARD};
static vector<bool> quadraticForms = {false, true};

BOOST_FIXTURE_TEST_SUITE(pixel_correction_suite, PixelFixture)

BOOST_DATA_TEST_CASE(
    pixel_correction_test,
    data::make(positions) * 
    data::xrange(1000.0, 20000.0, 1000.0) * 
    (data::make(iterations) ^ data::make(requiredAccuracies)) *
    data::make(numericMethods) *
    data::make(quadraticForms),
    position, cloudHeight, maxIterations, requiredAccuracy, numericMethod, quadraticForm)
{
  auto geos = geotransform.calcGeoCoordsFromPix(position.x + 0.5, position.y + 0.5);

  auto ellips = corrector.transformToEllipsCoordinates(geos);

  BOOST_TEST((!isnan(ellips.first) && !isnan(ellips.second)));

  auto cloudXYZ = calculateCloudPosition(ellips.second * M_PI / 180.0, ellips.first * M_PI / 180.0, cloudHeight, a, eSqr, central_m);
  auto simGeos = calculateGEOSCorrdsFromXYZ(cloudXYZ, a, satelliteHeight);

  BOOST_TEST((!isnan(simGeos.first) && !isnan(simGeos.second)));

  auto corrected = corrector.calculateNewCoordinates(simGeos, ellips, cloudHeight, requiredAccuracy, maxIterations);

  BOOST_TEST((!isnan(corrected.first) && !isnan(corrected.second)));

  auto correctedGeos = corrector.transformToGeosCoordinates(corrected);
  
  auto res = sqrt(pow(correctedGeos.first - geos.first, 2) + pow(correctedGeos.second - geos.second, 2)); //Indirect test

  BOOST_CHECK_MESSAGE(true, res);

  BOOST_TEST(res <= requiredAccuracy);
}

BOOST_AUTO_TEST_SUITE_END()