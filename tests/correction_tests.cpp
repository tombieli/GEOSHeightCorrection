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

static vector<GeoPosition> positions = {
    {18.4239, 33.9253}, // Cape Town
    {3.6947, 40.4177},  // Madrit
    {47.9142, 15.7839}, // Brasília
    {18.6453, 54.3475}, // Gdańsk
    {18.9333, 69.6667}  // Tromsø
};

static vector<int> iterations = {100, 1000};
static vector<double> requiredAccuracies = {10, 1};
static vector<NumericMethod> numericMethods = {NumericMethod::NETWON, NumericMethod::LEVENBERG_MARQUARD};
static vector<bool> quadraticForms = {false, true};



BOOST_FIXTURE_TEST_SUITE(correction_suite, Fixture)

BOOST_DATA_TEST_CASE(
    correction_test,
    data::make(positions) * 
    data::xrange(1000.0, 20000.0, 1000.0) * 
    (data::make(iterations) ^ data::make(requiredAccuracies)) *
    data::make(numericMethods) *
    data::make(quadraticForms),
    position, cloudHeight, maxIterations, requiredAccuracy, numericMethod, quadraticForm)
{
  auto orginalGeos = corrector.transformToGeosCoordinates(make_pair(position.lon, position.lat));

  BOOST_TEST((!isnan(orginalGeos.first) && !isnan(orginalGeos.second)));

  auto cloudXYZ = calculateCloudPosition(position.lat * M_PI / 180.0, position.lon * M_PI / 180.0, cloudHeight, a, eSqr, central_m);
  auto geos = calculateGEOSCorrdsFromXYZ(cloudXYZ, a, satelliteHeight);

  auto ellips = corrector.transformToEllipsCoordinates(geos);

  BOOST_TEST((!isnan(ellips.first) && !isnan(ellips.second)));

  auto corrected = corrector.calculateNewCoordinates(geos, ellips, cloudHeight, requiredAccuracy, maxIterations, numericMethod, quadraticForm);

  BOOST_TEST((!isnan(corrected.first) && !isnan(corrected.second)));

  corrected = corrector.transformToGeosCoordinates(corrected);

  auto res = sqrt(pow(orginalGeos.first - corrected.first, 2) + pow(orginalGeos.second - corrected.second, 2));

  BOOST_CHECK_MESSAGE(true, res);

  BOOST_TEST(res <= requiredAccuracy);
}

BOOST_AUTO_TEST_SUITE_END()