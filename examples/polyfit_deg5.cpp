#include <iostream>
#include <iomanip>
#include "tinyfit/polyfit.hpp"

// y = 1 + 2x − 0.5x² + 0.25x³ − 0.1x⁴ + 0.05x⁵  (exact)
constexpr std::array<double,6> xs{0, 1, 2, 3, 4, 5};
constexpr std::array<double,6> ys{
    1.0,
    2.6999999999999997,
    5.0,
    13.3,
    42.6,
    123.5
};

int main()
{
    constexpr auto c =
        tinyfit::polyfit::least_squares<5>(xs, ys);   // compile-time

    std::cout << std::setprecision(17);
    for (int i = 0; i < 6; ++i)
        std::cout << "a" << i << " = " << c[i] << '\n';
}