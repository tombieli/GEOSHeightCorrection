#include <iostream>
#include <string>
#include <ogr_spatialref.h>
#include <boost/program_options.hpp>

#include "program_options_ext.h"

namespace ksg {

    struct ProjSRS {
        OGRSpatialReference srs;
    };


    inline std::istream& operator>>(std::istream& in, ProjSRS& srs) {
        std::string line;
        std::getline(in, line);

        if(srs.srs.importFromProj4(line.c_str()) != OGRERR_NONE) {
            throw std::logic_error("Unable to parse proj srs: \""+line+"\"");
        }

        return in;
    }

    template<char SEP=','>
    struct Geotransform {
        double geotransform[6];

        std::pair<double,double> calcGeoCoordsFromPix(double xpix, double ypix)
        {
            double xgeo = geotransform[0] + geotransform[1] * xpix + geotransform[2] * ypix;
            double ygeo = geotransform[3] + geotransform[4] * xpix + geotransform[5] * ypix;
            return std::pair<double,double>(xgeo,ygeo);
        }
    };

    template<char SEP>
    inline std::istream& operator>>(std::istream& in, Geotransform<SEP>& geo) {
        std::tuple<double,double,double,double,double,double> packed;
        ksg::TupleParser::parse_tuple_with_separator<SEP, double,double,double,double,double,double>(in, packed);
        std::tie(geo.geotransform[0],geo.geotransform[1],geo.geotransform[2], 
                 geo.geotransform[3],geo.geotransform[4],geo.geotransform[5]) = packed;
        return in;
    }

    template<char SEP=','>
    struct ImageDimmensions {
        size_t x,y;
    };

    template<char SEP>
    inline std::istream& operator>>(std::istream& in, ImageDimmensions<SEP>& dims) {
        std::tuple<size_t,size_t> packed;
        ksg::TupleParser::parse_tuple_with_separator<SEP, size_t, size_t>(in, packed);
        std::tie(dims.x, dims.y) = packed;
        return in;
    }

}