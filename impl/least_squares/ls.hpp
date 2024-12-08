#ifndef LS_HPP
#define LS_HPP
#include "../math/pow_constexpr.hpp"
#include <array>
namespace tinyfit {
namespace impl {
namespace curve_fit {

// Computes the polynomial coefficients for least squares fitting
template <typename T, std::size_t degree, std::size_t data_size>
constexpr std::array<T, degree + 1>
fit(const std::array<T, data_size> &x_values,
    const std::array<T, data_size> &y_values) {
  std::array<T, 2 * degree + 1> X = {}; // Stores Σ(x^i) terms
  std::array<T, degree + 1> Y = {};     // Stores Σ(y * x^i) terms
  std::array<std::array<T, degree + 2>, degree + 1> B = {}; // Augmented matrix
  std::array<T, degree + 1> coefficients = {}; // Final coefficients

  // Calculate Σ(x^i) for i = 0 to 2*degree
  for (std::size_t i = 0; i < 2 * degree + 1; ++i) {
    for (std::size_t j = 0; j < data_size; ++j) {
      X[i] += math::pow_constexpr(x_values[j], static_cast<int>(i));
    }
  }

  // Build the normal matrix excluding the last column
  for (std::size_t i = 0; i <= degree; ++i) {
    for (std::size_t j = 0; j <= degree; ++j) {
      B[i][j] = X[i + j];
    }
  }

  // Calculate Σ(y * x^i) for i = 0 to degree
  for (std::size_t i = 0; i <= degree; ++i) {
    for (std::size_t j = 0; j < data_size; ++j) {
      Y[i] +=
          math::pow_constexpr(x_values[j], static_cast<int>(i)) * y_values[j];
    }
    B[i][degree + 1] = Y[i]; // Add Y as the last column of the augmented matrix
  }

  // Gaussian Elimination to solve the linear system
  std::size_t n = degree + 1;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = i + 1; k < n; ++k) {
      if (B[i][i] < B[k][i]) {
        for (std::size_t j = 0; j <= n; ++j) {
          std::swap(B[i][j], B[k][j]);
        }
      }
    }

    for (std::size_t k = i + 1; k < n; ++k) {
      T t = B[k][i] / B[i][i];
      for (std::size_t j = 0; j <= n; ++j) {
        B[k][j] -= t * B[i][j];
      }
    }
  }

  // Back substitution to calculate coefficients
  for (int i = n - 1; i >= 0; --i) {
    coefficients[i] = B[i][n];
    for (std::size_t j = i + 1; j < n; ++j) {
      coefficients[i] -= B[i][j] * coefficients[j];
    }
    coefficients[i] /= B[i][i];
  }

  return coefficients;
}
} // namespace curve_fit
} // namespace impl
} // namespace tinyfit
#endif