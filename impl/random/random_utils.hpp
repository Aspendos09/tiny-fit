#ifndef RANDOM_UTILS_HPP
#define RANDOM_UTILS_HPP
#include "CompileTimeRandom.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <type_traits>
namespace tinyfit {
namespace utils {
template <typename T, size_t size>
  requires std::is_arithmetic_v<T>
constexpr std::array<T, size> generate_random_numbers(T min, T max) {

  // this wrapper can be eliminated with modification of random number
  // generation library but. i dont wish to change implementation of random
  // number generation library right now. so this is not optimized wrapper but
  // at least it shouldn't affect application for compile time calculation.
  std::array<uint64_t, size> intermediate_vals{};

  Dynlec::CTRandomStream<size>::Call(
      [&intermediate_vals](uint64_t index, uint64_t n) {
        intermediate_vals[index] = n;
      });

  std::array<T, size> result{};

  for (auto idx{0}; idx < size; idx++) {

    const auto val = intermediate_vals[idx];

    long double normalized =
        static_cast<long double>(val) /
        static_cast<long double>(std::numeric_limits<uint64_t>::max());

    T scaled = static_cast<T>(min + normalized * (max - min));
    result[idx] = scaled;
  }

  return result;
}
} // namespace utils
} // namespace tinyfit

#endif