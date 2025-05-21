#include <iostream>
#include <iomanip>
#include "tinyfit/random.hpp"
#include "tinyfit/pso.hpp"
#include "tinyfit/polyfit.hpp"

// ----- raw data (6 points) -----
constexpr std::array<double,6> xs{0,1,2,3,4,5};
constexpr std::array<double,6> ys{1,2.7,5,13.3,42.6,123.5};

// ----- objective: RSS of degree-5 poly (6 coffs) -----
constexpr auto rss = [](const std::array<double,6>& c) {
    double err = 0.0;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        double x = xs[i];
        // Horner evaluation
        double yhat = c[5];
        for (int d = 4; d >= 0; --d) yhat = yhat * x + c[d];
        double diff = yhat - ys[i];
        err += diff * diff;
    }
    return err;
};

// ----- compile-time PSO (60×1000) -----
constexpr auto best =
    tinyfit::pso::run<double, 6, 60, 1000, rss>(
        tinyfit::random::constexpr_lcg{0xBEEF1234u});

int main() {
    std::cout << std::setprecision(17);
    std::cout << "compile_time_rss=" << best.second << '\n';
    for (int i = 0; i < 6; ++i)
        std::cout << "c" << i << '=' << best.first[i] << '\n';
}
