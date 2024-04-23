#pragma once
#include <iostream>
#include <string>
#include <string_view>

namespace ksg
{
    enum class NumericMethod {
        NETWON,
        LEVENBERG_MARQUARD
    };

    constexpr const char* NM_NETWON_NAME = "NETWON";
    constexpr const char* NM_NETWON_DESC = "Netwon numeric method. I could be lighter in computation effort. However it fails to perform correction in subsatellite point.";
    constexpr const char* NM_LEVENBERG_MARQUARD_NAME = "LEVENBERG_MARQUARD";
    constexpr const char* NM_LEVENBERG_MARQUARD_DESC = "Levenberg-Marquard numeric method. This method is more roboust. It allows for correction near egde of view disc.";

    inline static std::istream& operator>>(std::istream& in, NumericMethod& method) {
        std::string word;
        in >> word;

        if(word.compare(NM_NETWON_NAME) == 0) {
            method = NumericMethod::NETWON;
        } else if(word.compare(NM_LEVENBERG_MARQUARD_NAME)) {
            method = NumericMethod::LEVENBERG_MARQUARD;
        } else {
            throw std::logic_error("Unknown numeric method: \""+word+"\"");
        }
        return in;
    }

    inline static std::ostream& operator<<(std::ostream& out, const NumericMethod& method) {
        switch(method) {
            case NumericMethod::NETWON:
                out << NM_NETWON_NAME;
                break;
            case NumericMethod::LEVENBERG_MARQUARD:
                out << NM_LEVENBERG_MARQUARD_NAME;
                break;
            default:
                out << "Unknown method";
                break;
        }

        return out;
    }
}