#include <iostream>

namespace ksg
{
    struct GeoPosition
    {
        double lon, lat;
    };

    inline std::ostream &operator<<(std::ostream &out, const GeoPosition &pos)
    {
        out << "GeoPosition(lon: " << pos.lon << ", lat: " << pos.lat << ")";
        return out;
    }

    struct PixPosition
    {
        unsigned long x,y;
    };

    inline std::ostream &operator<<(std::ostream &out, const PixPosition &pos)
    {
        out << "PixPosition(x: " << pos.x << ", y: " << pos.y << ")";
        return out;
    }
}