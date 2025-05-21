#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace tinyfit::polyfit {

/* Simple Gauss–Jordan solver, O(N³) */
template <std::size_t N>
constexpr std::array<double, N> gauss_jordan(
    std::array<std::array<double, N + 1>, N> m) {

  for (std::size_t i = 0; i < N; ++i) {
    
    std::size_t piv = i;
    for (std::size_t k = i + 1; k < N; ++k)
      if (std::fabs(m[k][i]) > std::fabs(m[piv][i]))
        piv = k;
    std::swap(m[i], m[piv]);

    double div = m[i][i];
    for (double& v : m[i])
      v /= div;

    for (std::size_t r = 0; r < N; ++r) {
      if (r == i)
        continue;
      double f = m[r][i];
      for (std::size_t c = 0; c < N + 1; ++c)
        m[r][c] -= f * m[i][c];
    }
  }
  std::array<double, N> out{};
  for (std::size_t i = 0; i < N; ++i)
    out[i] = m[i][N];
  return out;
}

/* Ordinary least squares polynomial fit */
template <std::size_t Degree, std::size_t M>
constexpr std::array<double, Degree + 1> least_squares(
    const std::array<double, M>& xs, const std::array<double, M>& ys) {
  constexpr std::size_t N = Degree + 1;

  /* Σ xᵏ table */
  std::array<double, 2 * Degree + 1> X{};
  for (std::size_t k = 0; k < X.size(); ++k)
    for (double x : xs)
      X[k] += std::pow(x, k);

  /* Build augmented normal matrix */
  std::array<std::array<double, N + 1>, N> B{};
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j)
      B[i][j] = X[i + j];
    for (std::size_t m = 0; m < M; ++m)
      B[i][N] += std::pow(xs[m], i) * ys[m];
  }
  return gauss_jordan<N>(B);
}

}  // namespace tinyfit::polyfit
