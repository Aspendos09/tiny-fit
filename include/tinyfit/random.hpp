#pragma once
#include <cstdint>
#include <random>

namespace tinyfit::random {

/* ---------- constexpr 64-bit LCG ---------- */
struct constexpr_lcg {
  std::uint64_t state;

  constexpr explicit constexpr_lcg(std::uint64_t seed) : state(seed) {}

  static constexpr std::uint64_t next(std::uint64_t s) {
    return s * 6364136223846793005ull + 1ull;
  }

  constexpr double operator()() {
    state = next(state);
    return double((state >> 11) * 0x1p-53);  // [0,1)
  }
};

/* ---------- runtime MT19937-64 ---------- */
struct runtime_mt19937 {
  std::mt19937_64 eng{std::random_device{}()};
  std::uniform_real_distribution<double> dist{0.0, 1.0};

  double operator()() { return dist(eng); }
};

}  // namespace tinyfit::random
