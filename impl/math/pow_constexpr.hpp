#ifndef POW_CONSTEXPR_HPP
#define POW_CONSTEXPR_HPP
namespace tinyfit {
namespace math {

// Efficient, constexpr non-recursive power function
template <typename T> constexpr T pow_constexpr(T base, int exp) {
  T result = 1;
  while (exp > 0) {
    if (exp % 2 == 1) { // If exp is odd
      result *= base;
    }
    base *= base;
    exp /= 2;
  }
  return result;
}
} // namespace math
} // namespace tinyfit
#endif