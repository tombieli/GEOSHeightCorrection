#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>
#include "correction.h"
#include "georeference_utils.h"
#include "correction_utils.h"

using namespace std;

template <char SEP = ','>
struct Range
{
    size_t min, max;
};

template <char SEP>
inline std::istream &operator>>(std::istream &in, Range<SEP> &range)
{
    std::tuple<size_t, size_t> packed;
    ksg::TupleParser::parse_tuple_with_separator<SEP, size_t, size_t>(in, packed);
    std::tie(range.min, range.max) = packed;
    return in;
}

template <char SEP = ','>
struct HeightsRange
{
    double min, max, step;
};

template <char SEP>
inline std::istream &operator>>(std::istream &in, HeightsRange<SEP> &range)
{
    std::tuple<double, double, double> packed;
    ksg::TupleParser::parse_tuple_with_separator<SEP, double, double, double>(in, packed);
    std::tie(range.min, range.max, range.step) = packed;
    return in;
}

template <char SEP>
inline std::ostream &operator<<(std::ostream &out, const HeightsRange<SEP> &range)
{
    out << range.min << SEP << range.max << SEP << range.step;
    return out;
}

static const string nmDesc = string() +
        "Numeric method. Either:\n" +
        ksg::NM_NETWON_NAME + "\n" +
        ksg::NM_NETWON_DESC + "\n" +
        ksg::NM_LEVENBERG_MARQUARD_NAME + "\n" +
        ksg::NM_LEVENBERG_MARQUARD_DESC + "\n";

boost::program_options::variables_map init(int argc, char **argv)
{
    boost::program_options::options_description options("GEOS Height correction table generator options");
    options.add_options()
    ("help", "Shows help message")
    ("projection-description", boost::program_options::value<ksg::ProjSRS>(), "Description of geos projection in Proj format, used to generate tables.")
    ("geotransform", boost::program_options::value<ksg::Geotransform<>>(), "Six comma separated numbers describining raster location in SRS. See GDAL Raster Model.")
    ("image-dimensions", boost::program_options::value<ksg::ImageDimmensions<>>(), "Comma separated image dimmensions: width and height in pixels. See GDAL Raster Model.")
    ("columns-scope", boost::program_options::value<Range<>>(), "Comma separated, 1 based, inclusive scope for rasters columns. If specified, generation of table will be limited to those columns.")
    ("lines-scope", boost::program_options::value<Range<>>(), "Comma separated, 1 based, inclusive scope for rasters lines. If specified, generation of table will be limited to those lines.")
    ("heights-range", boost::program_options::value<HeightsRange<>>()->default_value({0.5, 20.0, 0.5}), "Range and of cloud heights. Should contain following comma separated values min, max, step. All numbers should be expressed in km.")
    ("output", boost::program_options::value<std::string>(), "Location of result table")
    ("requierd-accuracy", boost::program_options::value<double>()->default_value(10), "Required accuracy [m].")
    (
        "numeric-method", 
        boost::program_options::value<ksg::NumericMethod>()->default_value(ksg::NumericMethod::LEVENBERG_MARQUARD),
        nmDesc.c_str()
    )
    ("use-squared-target", "Use squared targed function.")
    ("iterations-limit", boost::program_options::value<int>()->default_value(100), "Maximum number of iteration per pixel.");

    auto ret = boost::program_options::variables_map();

    auto style = boost::program_options::command_line_style::unix_style;
    style = (decltype(style))(style ^ boost::program_options::command_line_style::allow_short);
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options, style), ret);
    boost::program_options::notify(ret);

    if (argc <= 1 || ret.count("help") >= 1)
    {
        std::cout << options << std::endl;
        exit(1);
    }

    return ret;
}

void check_required_option(const boost::program_options::variables_map &variablesMap, const std::string &key)
{
    if (variablesMap.count(key) == 0)
    {
        throw logic_error("Please provide --" + key);
    }
}

int main(int argc, char **argv)
{
    auto variablesMap = init(argc, argv);

    check_required_option(variablesMap, "projection-description");
    auto srs = variablesMap["projection-description"].as<ksg::ProjSRS>();

    auto projectionName = string(srs.srs.GetAttrValue("PROJECTION"));

    if (projectionName != "Geostationary_Satellite")
        throw runtime_error("Height correction requires Geostationary_Satellite projection.");

    check_required_option(variablesMap, "geotransform");
    auto geo = variablesMap["geotransform"].as<ksg::Geotransform<>>();

    check_required_option(variablesMap, "image-dimensions");
    auto dims = variablesMap["image-dimensions"].as<ksg::ImageDimmensions<>>();

    size_t startX, endX;
    size_t startY, endY;

    if (variablesMap.count("columns-scope") > 0)
    {
        auto customScope = variablesMap["columns-scope"].as<Range<>>();
        startX = customScope.min - 1;
        endX = customScope.max - 1;
    }
    else
    {
        startX = 0;
        endX = dims.x - 1;
    }

    if (variablesMap.count("lines-scope") > 0)
    {
        auto customScope = variablesMap["lines-scope"].as<Range<>>();
        startY = customScope.min - 1;
        endY = customScope.max - 1;
    }
    else
    {
        startY = 0;
        endY = dims.y - 1;
    }

    auto heightRange = variablesMap["heights-range"].as<HeightsRange<>>();

    auto requiredAccuracy = variablesMap["requierd-accuracy"].as<double>();

    auto iterationLimit = variablesMap["iterations-limit"].as<int>();

    auto numericMethod = variablesMap["numeric-method"].as<ksg::NumericMethod>();

    auto useQuadraticForm = variablesMap.count("use-squared-target") > 0;

    check_required_option(variablesMap, "output");
    auto outputName = variablesMap["output"].as<std::string>();

    std::ofstream output(outputName);

    output.exceptions(std::ios::failbit | std::ios::badbit);

    GEOSHeightCorrector corrector(srs.srs);

    std::vector<double> heights;

    for(double h = heightRange.min; h <= heightRange.max; h += heightRange.step) {
        heights.push_back(h * 1000);    
    }

    for (size_t y = startY; y <= endY; y++)
    {
        for (size_t x = startX; x <= endX; x++)
        {
            auto geosCoords = geo.calcGeoCoordsFromPix(x + 0.5, y + 0.5);
            auto ellipsCoords = corrector.transformToEllipsCoordinates(geosCoords);

            if(isnan(ellipsCoords.first) || isnan(ellipsCoords.second)) {
                std::cerr << "Pixel " << x + 1 <<", " << y + 1 << " is out of scope - skipping.\n";
                continue;
            }
            
            std::vector<std::pair<double,double>> results;
            results.resize(heights.size());

            #pragma omp parallel for shared(results) 
            for(size_t hi = 0; hi < heights.size(); hi++) {
                auto h = heights[hi];
                results[hi] = corrector.calculateNewCoordinates(geosCoords, ellipsCoords, h, requiredAccuracy, iterationLimit, numericMethod, useQuadraticForm);
            }

            if(std::all_of(results.begin(), results.end(), [](auto pair)  { return isnan(pair.first) && isnan(pair.second); } )) {
                std::cerr << "Failed to generate table for pixel " << x + 1 <<", " << y + 1 << ".\n";
                continue;
            }

            output << y + 1 << " " << x + 1;

            for(size_t hi = 0; hi < results.size(); hi++) {
                auto res = results[hi];
                output << " " << res.second << " " << res.first;
            }

            output << "\n";
        }
    }

    return 0;
}
