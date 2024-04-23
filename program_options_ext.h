#include <tuple>
#include <iostream>

namespace ksg
{

    class TupleParser
    {
        template <int n, int rest, char SEP, typename... T>
        class Helper
        {
            public:
            static inline void parse_tuple_with_separator_(std::istream &in, std::tuple<T...> &tuple)
            {
                in >> std::get<n>(tuple);
                char sep;
                in >> sep;
                if (sep != SEP)
                {
                    char expected = SEP;
                    throw std::logic_error("Bad separator, expected \"" + std::string((const char *)&expected, 1) + "\", but was \"" + std::string((const char *)&sep, 1) + "\"");
                }
                Helper<n + 1, rest - 1, SEP, T...>::parse_tuple_with_separator_(in, tuple);
            }
        };

        template<int n, char SEP, typename... T>
        class Helper<n, 1, SEP, T...>
        {
            public:
            static inline void parse_tuple_with_separator_(std::istream &in, std::tuple<T...> &tuple) {
                in >> std::get<n>(tuple);
            }
        };

    public:
        template <char SEP, typename... T>
        static inline void parse_tuple_with_separator(std::istream &in, std::tuple<T...> &tuple)
        {
            Helper<0, sizeof...(T), SEP, T...>::parse_tuple_with_separator_(in, tuple);
        }
    };

}