#include <cmath>

namespace ksg
{
    constexpr static inline double N(double phi_earth, double a, double eSqr)
    {
        using namespace std;
        return a / sqrt(1 - eSqr * pow(sin(phi_earth), 2));
    }

    constexpr static inline double dNdphi(double phi_earth, double a, double eSqr)
    {
        using namespace std;
        return (a * eSqr * sin(phi_earth) * cos(phi_earth)) / pow(1 - eSqr * pow(sin(phi_earth), 2), 3 / 2);
    }

}